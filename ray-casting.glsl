
const int MAX_MARCHING_STEPS = 10;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
const int SAMPLES = 16;
const float PI = 3.1415926535;


vec3 circle_intersect(vec2 o, vec2 d, vec2 c, float r)
{
	vec2 u = c-o;
	float l = dot(u,d);
	float d2 = dot(u,u ) - l*l;
  float k = r*r - d2;
	if (k<0.0)
		return vec3(-1.0,0.0,0.0);
	float s = sqrt(k);
	float t1 = l - s;
  float t2 = l + s;
	float t3 = t1>0.0?t1:t2;
  vec2 n = ((o+t3*d) - c)/r;
  return vec3(t3,n);
}

vec2 light;

void move_light()
{
  float a = 3.1415926 * iTime / 5.0;
  light = vec2(cos(a),sin(a))*0.7 ;//+ vec2(0.3,-0.5);
}


vec4 scene_intersect(vec2 o, vec2 d) {
  vec3 s = vec3(1e6,0.0,0.0);
  float li = -1.0;
  vec3 i = circle_intersect(o,d,vec2(0.,0.),0.2);
  if(i.x>EPSILON&&i.x<s.x) {s=i;li=0.0;}
  i = circle_intersect(o,d,light,0.3);
  if(i.x>EPSILON&&i.x<s.x) {s=i;li=1.0;}
  i = circle_intersect(o,d,-light,0.3);
  if(i.x>EPSILON&&i.x<s.x) {s=i;li=0.0;}
  return vec4(s,li);
}



bool is_total_reflection(vec2 i,vec2 n,float refrac){
  float cosp = dot(i,n);
  float cosp22 = 1.0 - refrac*refrac*(1.0 - cosp*cosp);
  return cosp22<0.0;
}

vec2 refractR(vec2 i,vec2 n,float refrac){
  float cosp = dot(i,n);
  float cosp22 = 1.0 - refrac*refrac*(1.0 - cosp*cosp);
  float sc = cosp*refrac - sqrt(cosp22);
  return i * refrac + sc * n;
}

vec2 reflectR(vec2 i,vec2 n)
{
  vec2 ni = n*dot(n,i)*2.0;
  return i-ni;
}

vec3 ray_cast(in vec2 org,in vec2 odir) {
    float t = 0.0;
    vec2 d = odir;
    vec3 total = vec3(0.0);
    vec2 o = org;
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        vec4 r = scene_intersect(o,d);
        if(r.w < 0.0) return total;
        total += vec3(r.w);
        //return total;
        if(r.w > EPSILON) return total;
        o = o+d*r.x;
        vec2 n = r.yz;
        float cosq = dot(n,d);
        float refrac = cosq > 0.0 ? 1.3 : 1.0/1.3;
        n = cosq > 0.0 ? -n : n;
        d = refract(d, n, refrac);

        if(dot(d,d)<0.0)
           d  = reflect(d, n);
        o += d*EPSILON*100.0;

    }
    return total;
}


float hash1(inout float seed) {
    return fract(sin(seed += 0.1)*43758.5453123);
}

vec2 hash2(inout float seed) {
  return fract(sin(vec2(seed+=0.1,seed+=0.1))*vec2(43758.5453123,22578.1459123));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 q = fragCoord.xy / iResolution.xy;
    vec2 p = -1.0 + 2.0 * (fragCoord.xy) / iResolution.xy;
    p.x *= iResolution.x/iResolution.y;
    float seed = p.x + p.y * 3.43121412313 + fract(1.12345314312*iTime);
    vec3 total = vec3(0.0);
    move_light();
    for(int i=0;i<SAMPLES;i++)
    {
        float angle = (hash1(seed) + float(i))/float(SAMPLES)*PI*2.0 ;

        vec2 dir = vec2(cos(angle),sin(angle));
        total += ray_cast(p, dir);
        seed = mod( seed*1.1234567893490423, 13. );
    }
    total/=float(SAMPLES);
    total = pow(clamp(total,0.0,1.0),vec3(0.45));
    fragColor = vec4(total, 1.0);
}
