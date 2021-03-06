#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <vector>
#include <memory>
#include <cstdio>

#define SVPNG_OUTPUT std::fstream& os
#define SVPNG_PUT(u) os.put(u)
#include "svpng\svpng.h"
#undef SVPNG_OUTPUT
#undef SVPNG_PUT



constexpr int nthread = 5;
using unit = float;
constexpr unit TWO_PI = 6.28318530718;
constexpr unit PI = 3.1415926535;
constexpr size_t sample_count = 32;
using BYTE = unsigned char;
constexpr size_t width = 1080, height = 1080;
struct emissive
{
	BYTE r, g, b;
}img[width][height], DF_img[width][height];

struct urgb
{
	unit r, g, b;
	urgb& operator+=(const urgb& o)
	{
		r += o.r;
		g += o.g;
		b += o.b;
		return *this;
	}
	urgb& operator*=(const urgb& o)
	{
		r *= o.r;
		g *= o.g;
		b *= o.b;
		return *this;
	}
	urgb& operator*=(unit o)
	{
		r *= o;
		g *= o;
		b *= o;
		return *this;
	}
};



unit circleSDF(unit x, unit y, unit cx, unit cy, unit r) {
	unit ux = x - cx, uy = y - cy;
	return sqrt(ux * ux + uy * uy) - r;
}


struct hit_record
{
	bool hit;
	bool in;
	unit t;
	unit x, y;
	unit dx, dy;
};

hit_record circle_intersect(unit ox, unit oy, unit dx, unit dy, unit cx, unit cy, unit r)
{
	hit_record rec;
	unit ux = cx - ox;
	unit uy = cy - oy;
	unit us = ux*ux + uy*uy;
	unit l = ux*dx + uy*dy;
	unit d2 = us - l*l;
	if (d2 > r*r)
	{
		rec.hit = false;
		return rec;
	}
	unit s = std::sqrt(r*r - d2);
	unit t = l - s;
	rec.in = false;
	if (t < 0)
	{
		t = l + s;
		rec.in = true;
	}
	if (t < 0)
	{
		rec.hit = false;
		return rec;
	}
	rec.t = t;
	rec.x = ox + dx*t;
	rec.y = oy + dy*t;
	rec.dx = (rec.x - cx) / r;
	rec.dy = (rec.y - cy) / r;
	rec.hit = true;
	return rec;
}

struct material
{
	urgb basecolor;
	urgb emissive = { 0,0,0 };
	float reflectivity = 1;
	float transparency = 1;
	float refractive = 1.3;
};

struct geo
{
	material _material;
	virtual unit SDF(unit x, unit y) = 0;
	virtual hit_record intersect(unit ox, unit oy, unit dx, unit dy) = 0;
	geo(material e) :_material(e) {}
};

struct circle : public geo
{
	unit _cx, _cy, _r;
	virtual unit SDF(unit x, unit y) override { return circleSDF(x, y, _cx, _cy, _r); }
	virtual hit_record intersect(unit ox, unit oy, unit dx, unit dy) override { return circle_intersect(ox, oy, dx, dy, _cx, _cy, _r); }
	circle(unit x, unit y, unit r, material e) :_cx(x), _cy(y), _r(r), geo(e) {}
};

std::vector<material> materials;
std::vector<std::shared_ptr<geo>> scene;

unit static_DF[width][height];
int geo_id[width][height];



void build_DF()
{
	scene.push_back(std::make_shared<circle>(0.3, 0.3, 0.10, material{ urgb{ 0.f, 0.f, 0.f } }));
	scene.push_back(std::make_shared<circle>(0.5, 0.5, 0.07, material{ urgb{ 2., 2., 2. }, urgb{ 2., 2., 2. }, 0., 0., 0. }));
	scene.push_back(std::make_shared<circle>(0.3, 0.7, 0.15, material{ urgb{ 0.2f, 0.5f, 0.8f } }));
	scene.push_back(std::make_shared<circle>(0.8, 0.7, 0.07, material{ urgb{ 0.5f, 0.1f, 0.8f } }));
	scene.push_back(std::make_shared<circle>(0.5, 0.3, 0.06, material{ urgb{ 0.7f, 0.5f, 0.1f } }));
	scene.push_back(std::make_shared<circle>(0.7, 0.6, 0.05, material{ urgb{ 0.2f, 0.8f, 0.3f } }));
	scene.push_back(std::make_shared<circle>(0.8, 0.5, 0.06, material{ urgb{ 0.7f, 0.1f, 0.3f } }));
	//scene.push_back(std::make_shared<circle>(0.8, 0.9, 0.05, material{ urgb{ 1.f, 1.2f, 2.5f }, urgb{ 1.f, 1.2f, 2.5f }, 0, 0, 0 }));
	scene.push_back(std::make_shared<circle>(0.8, 0.25, 0.10, material{ urgb{ 0.7f, 0.7f, 0.7f } }));

	//scene.push_back(std::make_shared<circle>(0.5, -100, 100.1, material{ urgb{ 0.1f, 0.2f, 0.35f } }));
	int size = scene.size();
	
	std::fstream os;
	os.open("DF.png", std::fstream::out | std::fstream::binary);

	materials.resize(size);
	for (int i = 0; i < size; i++)
	{
		materials[i] = scene[i]->_material;
	}
#pragma omp parallel for schedule(dynamic, 1) /*private(sample)*/ num_threads(nthread)
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			unit sd = 10.;
			int id = -1;
			for (int i = 0; i < size; i++)
			{
				unit s = scene[i]->SDF(static_cast<unit>(x) / width, static_cast<unit>(y) / height);
				if (id == -1 || s < sd)
				{
					sd = s;
					id = i;
				}
			}

			static_DF[x][y] = sd;
			geo_id[x][y] = id;

			auto clamp = [&](unit v) {return std::fmin(std::fmax(v, 0.)*255., 255.); };
			auto& pixel = DF_img[x][y];
			pixel.r = static_cast<BYTE>(clamp(sd));
			pixel.g = static_cast<BYTE>(clamp(sd));
			pixel.b = static_cast<BYTE>(clamp(sd));
		}
	}

	svpng(os, width, height, reinterpret_cast<const BYTE*>(DF_img), 0);
}

constexpr unit pixel_w = 1. / width;
constexpr unit pixel_h = 1. / height;

struct scene_pixel
{
	unit sd;
	material mat;
};


auto sample_static_DF(unit x, unit y) {
	

	auto clamp = [](unit v) -> unit
	{
		return std::fmin(std::fmax(v, 0), 0.9999);
	};

	int ix = static_cast<int>(clamp(x)*width);
	int iy = static_cast<int>(clamp(y)*height);
	

	return scene_pixel{ static_DF[ix][iy],materials[geo_id[ix][iy]] };
}

auto sample_dynamic_DF(unit x, unit y) {
	int size = scene.size();
	unit sd = 10.;
	int id = -1;
	for (int i = 0; i < size; i++)
	{
		unit s = scene[i]->SDF(x, y);
		if (id == -1 || s < sd)
		{
			sd = s;
			id = i;
		}
	}
	return scene_pixel{ sd, scene[id]->_material };
}

constexpr unit static_DF_EPSILON = pixel_w < pixel_h ? pixel_w : pixel_h;
constexpr unit dynamic_DF_EPSILON = 1e-6;

constexpr auto sample_DF = sample_dynamic_DF;
constexpr unit DF_EPSILON = sample_DF == sample_dynamic_DF ? dynamic_DF_EPSILON : static_DF_EPSILON;

inline auto sample_DF_vector(unit x,unit y)
{
	struct
	{
		unit x, y;
	}vec;
	vec.x = (sample_DF(x + DF_EPSILON, y).sd - sample_DF(x - DF_EPSILON, y).sd)*(0.5 / DF_EPSILON);
	vec.y = (sample_DF(x, y + DF_EPSILON).sd - sample_DF(x, y - DF_EPSILON).sd)*(0.5 / DF_EPSILON);
	return vec;
}

auto reflect(unit ix, unit iy, unit nx, unit ny) {
	unit idotn2 = (ix * nx + iy * ny) * 2.0f;
	struct {
		unit x, y;
	}vec;
	vec.x = ix - idotn2 * nx;
	vec.y = iy - idotn2 * ny;
	return vec;
}

auto is_total_reflection(unit ix, unit iy, unit nx, unit ny, unit refrac) {
	unit cosp = ix*nx + iy*ny;
	unit cosp22 = 1 - refrac*refrac*(1 - cosp*cosp);
	unit sc = cosp*refrac - sqrt(cosp22);
	return cosp22 < 0;
}


auto refract(unit ix,unit iy,unit nx,unit ny,unit refrac){
	struct {
		unit x, y;
	}vec;
	unit cosp = ix*nx + iy*ny;
	unit cosp22 = 1 - refrac*refrac*(1 - cosp*cosp);
	unit sc = cosp*refrac - sqrt(cosp22);
	vec.x = ix*refrac + sc * nx;
	vec.y = iy*refrac + sc * ny;
	return vec;
}

urgb trace_ray_marching(unit ox, unit oy, unit a, int depth = 0) {
    constexpr size_t MAX_depth = 5;
	constexpr size_t MAX_STEP = 40;
	constexpr unit MAX_DISTANCE = 5.0;
	if (depth > MAX_depth) return { 0,0,0 };
	unit dx = std::cos(a), dy = std::sin(a);
	unit edge = DF_EPSILON;
	auto r = sample_DF(ox, oy);
	r.sd = std::fabs(r.sd);

	if (depth == 0 && r.sd < DF_EPSILON) {
		if(r.mat.transparency == 0) return r.mat.emissive;
	}

	

	unit escape = DF_EPSILON * 10;
	unit t = 0.;
	for (int i = 0; i < MAX_STEP && t < MAX_DISTANCE; i++) {
		if (r.sd < DF_EPSILON) {
			auto vec = sample_DF_vector(ox, oy);
			auto refrac = r.mat.refractive;
			auto cosq = vec.x*dx + vec.y*dy;
			if (cosq > 0)
			{
				refrac = 1 / refrac;
				vec.x = -vec.x;
				vec.y = -vec.y;
			}
			if (r.mat.transparency > 0 && !is_total_reflection(dx, dy, vec.x, vec.y, refrac))
			{
				auto out = refract(dx, dy, vec.x, vec.y, refrac);
				ox += out.x * escape;
				oy += out.y * escape;
				a = std::atan2(out.y, out.x);
				auto color = trace_ray_marching(ox, oy, a, depth + 1);
				if(cosq>0) color *= r.mat.basecolor;
				color += r.mat.emissive;
				return color;
			}
			else
			{
				if (std::rand() / RAND_MAX < r.mat.reflectivity)
				{
					auto out = reflect(dx, dy, vec.x, vec.y);
					ox += out.x * escape;
					oy += out.y * escape;
					a = std::atan2(out.y, out.x);
					auto color = trace_ray_marching(ox, oy, a, depth + 1);
					color += r.mat.emissive;
					return color;
				}
				else
				{
					a = std::atan2(vec.y, vec.x);
					a = a - PI / 2 + PI * static_cast<unit>(rand()) / RAND_MAX;
					dx = std::cos(a), dy = std::sin(a);
					ox += dx * escape;
					oy += dy * escape;
					auto color = trace_ray_marching(ox, oy, a, depth + 1);
					color *= r.mat.basecolor;
					color += r.mat.emissive;
					return color;
				}
			}
			
		}
		else {
			ox += dx * r.sd;
			oy += dy * r.sd;
			t += r.sd;
		}
		r = sample_DF(ox, oy);
		r.sd = std::fabs(r.sd);
	}
	return { 0,0,0 };
}





struct scene_record
{
	hit_record hit;
	material mat;
};

auto trace(unit ox, unit oy, unit dx, unit dy)
{
	int size = scene.size();
	hit_record rec;
	scene_record r;
	rec.hit = false;
	rec.t = 10.;
	r.hit = rec;
	int id = -1;
	for (int i = 0; i < size; i++)
	{
		hit_record r = scene[i]->intersect(ox,oy,dx,dy);
		if (r.hit && r.t < rec.t)
		{
			rec = r;
			id = i;
		}
	}
	if (id == -1)
	{
		return r;
	}
	r.hit = rec;
	r.mat = materials[id];
	return r;
}

urgb trace_ray_casting(unit ox, unit oy, unit a, int depth = 0) {

	constexpr size_t MAX_depth = 10;
	constexpr unit EPSILON = 1e-6;
	constexpr unit escape = EPSILON * 10;
	if (depth > MAX_depth) return { 0,0,0 };
	unit dx = std::cos(a), dy = std::sin(a);
	auto rec = trace(ox, oy, dx, dy);
	if (!rec.hit.hit) return { 0,0,0 };
	if (rec.hit.in)
	{
		if (rec.mat.transparency == 0) return rec.mat.emissive;
	}
	ox = rec.hit.x;
	oy = rec.hit.y;

	auto refrac = rec.mat.refractive;
	auto cosq = rec.hit.dx*dx + rec.hit.dy*dy;
	if (cosq > 0)
	{
		refrac = 1 / refrac;
		rec.hit.dx = -rec.hit.dx;
		rec.hit.dy = -rec.hit.dy;
	}
	if (rec.mat.transparency > 0 && !is_total_reflection(dx, dy, rec.hit.dx, rec.hit.dy, refrac))
	{
		auto out = refract(dx, dy, rec.hit.dx, rec.hit.dy, refrac);
		ox += out.x * escape;
		oy += out.y * escape;
		a = std::atan2(out.y, out.x);
		auto color = trace_ray_marching(ox, oy, a, depth + 1);
		if (cosq>0) color *= rec.mat.basecolor;
		color += rec.mat.emissive;
		return color;
	}
	else
	{
		if (std::rand() / RAND_MAX < rec.mat.reflectivity)
		{
			auto out = reflect(dx, dy, rec.hit.dx, rec.hit.dy);
			ox += out.x * escape;
			oy += out.y * escape;
			a = std::atan2(out.y, out.x);
			auto color = trace_ray_casting(ox, oy, a, depth + 1);
			color += rec.mat.emissive;
			return color;
		}
		else
		{
			a = std::atan2(rec.hit.dy, rec.hit.dx);
			a = a - PI / 2 + PI * static_cast<unit>(rand()) / RAND_MAX;
			dx = std::cos(a), dy = std::sin(a);
			ox += dx * escape;
			oy += dy * escape;
			auto color = trace_ray_casting(ox, oy, a, depth + 1);
			color *= rec.mat.basecolor;
			color += rec.mat.emissive;
			return color;
		}
	}
	
}



int main()
{
	std::fstream os;
	os.open("result.png", std::fstream::out | std::fstream::binary);
	urgb sample;
	std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
	build_DF();
	std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
	std::cout << "构建距离场 耗时: "
		<< ::std::chrono::duration_cast<::std::chrono::milliseconds>(
			end - start)
		.count()
		<< "\n";
	start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for schedule(dynamic, 1) private(sample) num_threads(nthread)
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			
			sample = { 0., 0., 0. };
			for (int s4 = 0; s4 < 4; s4++)
			{
				int sx = s4 % 2, sy = s4 / 2;
				unit r1 = 2. * (unit)rand() / RAND_MAX, dx = r1 < 1. ? sqrt(r1) - 1. : 1. - sqrt(2. - r1);
				unit r2 = 2. * (unit)rand() / RAND_MAX, dy = r2 < 1. ? sqrt(r2) - 1. : 1. - sqrt(2. - r2);
				unit xx = ((sx + .5 + dx) / 2 + x) / width;
				unit yy = ((sy + .5 + dy) / 2 + y) / height;
				for (int i = 0; i < sample_count; i++) {
					unit a = TWO_PI * (i + (unit)rand() / RAND_MAX) / sample_count;
					sample += trace_ray_marching(xx, yy, a);
				}
			}
			auto clamp = [&](unit v) {return std::fmin((v/(sample_count*4))*255., 255.); };
			auto& pixel = img[x][y];
			pixel.r = static_cast<BYTE>(clamp(sample.r));
			pixel.g = static_cast<BYTE>(clamp(sample.g));
			pixel.b = static_cast<BYTE>(clamp(sample.b));
		}
	}
	end = std::chrono::high_resolution_clock::now();
	std::cout << "渲染 耗时: "
		<< ::std::chrono::duration_cast<::std::chrono::milliseconds>(
			end - start)
		.count()
		<< "\n";

	svpng(os, width, height, reinterpret_cast<const BYTE*>(img), 0);
    return 0;
}

