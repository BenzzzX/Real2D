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



constexpr int nthread = 6;
using unit = float;
constexpr unit TWO_PI = 6.28318530718;
constexpr unit PI = 3.1415926535;
constexpr size_t sample_count = 2048;
using BYTE = unsigned char;
constexpr size_t width = 1080, height = 1080;
struct rgb
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
	rec.hit = true;
	return rec;
}

struct geo
{
	urgb _emissive;
	virtual unit SDF(unit x, unit y) = 0;
	virtual hit_record intersect(unit ox, unit oy, unit dx, unit dy) = 0;
	geo(urgb e) :_emissive(e) {}
};

struct circle : public geo
{
	unit _cx, _cy, _r;
	virtual unit SDF(unit x, unit y) override { return circleSDF(x, y, _cx, _cy, _r); }
	virtual hit_record intersect(unit ox, unit oy, unit dx, unit dy) override { return circle_intersect(ox, oy, dx, dy, _cx, _cy, _r); }
	circle(unit x, unit y, unit r, urgb e) :_cx(x), _cy(y), _r(r), geo(e) {}
};

std::vector<urgb> emissive;
std::vector<std::shared_ptr<geo>> scene;

unit static_DF[width][height];
int geo_id[width][height];



void build_DF()
{
	scene.push_back(std::make_shared<circle>(0.3, 0.3, 0.10, urgb{ 2., 2., 2. })); 
	scene.push_back(std::make_shared<circle>(0.5, 0.5, 0.07, urgb{ 0.f, 0.f, 0.f }));
	scene.push_back(std::make_shared<circle>(0.3, 0.7, 0.05, urgb{ 0.1f, 0.25f, 0.4f }));
	scene.push_back(std::make_shared<circle>(0.5, 0.7, 0.07, urgb{ 0.25f, 0.05f, 0.4f }));
	scene.push_back(std::make_shared<circle>(0.5, 0.3, 0.06, urgb{ 0.35f, 0.25f, 0.05f }));
	scene.push_back(std::make_shared<circle>(0.7, 0.6, 0.05, urgb{ 0.1f, 0.4f, 0.15f }));
	scene.push_back(std::make_shared<circle>(0.8, 0.5, 0.06, urgb{ 0.35f, 0.05f, 0.15f }));
	scene.push_back(std::make_shared<circle>(0.8, 0.9, 0.05, urgb{ 1.f, 1.2f, 2.5f }));
	scene.push_back(std::make_shared<circle>(0.9, 0.2, 0.08, urgb{ 0.1f, 0.2f, 0.35f }));
	int size = scene.size();
	
	std::fstream os;
	os.open("DF.png", std::fstream::out | std::fstream::binary);

	emissive.resize(size);
	for (int i = 0; i < size; i++)
	{
		emissive[i] = scene[i]->_emissive;
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
	urgb emissive;
};


auto sample_static_DF(unit x, unit y) {
	

	auto clamp = [](unit v) -> unit
	{
		return std::fmin(std::fmax(v, 0), 0.9999);
	};

	int ix = static_cast<int>(clamp(x)*width);
	int iy = static_cast<int>(clamp(y)*height);
	

	return scene_pixel{ static_DF[ix][iy],emissive[geo_id[ix][iy]] };
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
	return scene_pixel{ sd, scene[id]->_emissive };
}

auto sample_dynamic_DF_vector(unit x, unit y)
{
	unit EPSILON = 1e-6;
	unit dx = sample_dynamic_DF(x + EPSILON, y).sd - sample_dynamic_DF(x, y).sd;
	unit dy = sample_dynamic_DF(x, y + EPSILON).sd - sample_dynamic_DF(x, y).sd;
	unit l = std::sqrt(dx*dx + dy*dy);
	struct {
		unit x, y;
	} r{ dx / l, dy / l };
	return r;
}

urgb trace_ray_marching(unit ox, unit oy, unit a) {
	
	unit dx = std::cos(a), dy = std::sin(a);
	constexpr size_t MAX_STEP = 35;
	constexpr unit MAX_DISTANCE = 5.0;
	unit EPSILON = std::fmin(pixel_w, pixel_h);
	unit edge = EPSILON*3;
	auto r = sample_dynamic_DF(ox, oy);
	if (r.sd < EPSILON) {
#if 1
		return r.emissive;
#else
		if(r.sd < EPSILON - edge) return r.emissive;
		auto vec = sample_dynamic_DF_vector(ox, oy);
		auto dot = dx*vec.x + dy*vec.y;
		//c^2=a^2+b^2-2abcosC
		if (dot < 0) return r.emissive;
		do
		{
			auto step = std::fmax(std::abs(r.sd), EPSILON);
			ox += dx * step;
			oy += dy * step;
			r = sample_dynamic_DF(ox, oy);
		} while (r.sd < EPSILON); //escape
#endif
	}
	unit escape = EPSILON * 10;
	unit t = 0.;
	urgb acc = {0., 0., 0.};
	for (int i = 0; i < MAX_STEP && t < MAX_DISTANCE; i++) {
		if (r.sd < EPSILON) {
			a = a + PI/2 + PI * static_cast<unit>(rand()) / RAND_MAX;
			dx = std::cos(a), dy = std::sin(a);
			r.emissive *= 0.4;
			acc += r.emissive;
			ox += dx * escape;
			oy += dy * escape;
			t += escape;
		}
		else {
			ox += dx * r.sd;
			oy += dy * r.sd;
			t += r.sd;
		}
		r = sample_static_DF(ox, oy);
	}
	return acc;
}





struct scene_record
{
	unit x, y;
	bool hit;
	bool in;
	urgb emissive;
};

scene_record trace(unit ox, unit oy, unit dx, unit dy)
{
	int size = scene.size();
	hit_record rec;
	rec.hit = false;
	rec.t = 10.;
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
		return { rec.x,rec.y,false,false,{} };
	}
	else return { rec.x,rec.y,true,rec.in,emissive[id] };
}

urgb trace_ray_casting(unit ox, unit oy, unit a) {
	unit dx = std::cos(a), dy = std::sin(a);
	urgb acc = { 0., 0., 0. };
	constexpr size_t MAX_STEP = 10;
	unit EPSILON = std::fmin(pixel_w, pixel_h) * 2;
	for (int i = 0; i < MAX_STEP; i++) {
		auto rec = trace(ox, oy, dx, dy);
		if (!rec.hit) break;
		if (rec.in)
		{
			acc += rec.emissive;
			break;
		}
		rec.emissive *= 0.4;
		acc += rec.emissive;
		
		a = a + PI/2 + PI * static_cast<unit>(rand()) / RAND_MAX;
		dx = std::cos(a), dy = std::sin(a);
		ox = rec.x; oy = rec.y;
		ox += dx * EPSILON;
		oy += dy * EPSILON;
	}
	return acc;
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
					sample += trace_ray_casting(xx, yy, a);
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

