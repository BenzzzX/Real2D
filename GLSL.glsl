
const int MAX_MARCHING_STEPS = 100;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
const int SAMPLES = 32;
const float PI = 3.1415926535;

float sphereSDF(vec2 samplePoint,vec2 c,float r) {
    return length(samplePoint - c) - r;
}

vec2 unionR(vec2 a,vec2 b)
{
    return a.x<b.x?a:b;
}

vec2 intersectOp(vec2 a, vec2 b) {
    vec2 r = a.x > b.x ? b : a;
    r.x = a.x > b.x ? a.x : b.x;
    return r;
}

vec2 subtractOp(vec2 a, vec2 b) {
    vec2 r = a;
    r.x = (a.x > -b.x) ? a.x : -b.x;
    return r;
}

vec2 light;

void move_light()
{
  float a = 3.1415926 * iTime / 5.0;
  light = vec2(cos(a),sin(a))*0.9 ;//+ vec2(0.3,-0.5);
}


vec2 sceneSDF(vec2 samplePoint) {
  vec2 s1 = vec2(sphereSDF(samplePoint,light,0.25),1.0);
  vec2 s2 = vec2(sphereSDF(samplePoint,vec2(0.6,0.4),0.2),0.0);
  vec2 s3 = vec2(sphereSDF(samplePoint,vec2(0.6,0.4),0.11),0.0);
  vec2 s4 = vec2(sphereSDF(samplePoint,vec2(-0.77,-0.),0.4),0.0);
  vec2 s5 = vec2(sphereSDF(samplePoint,vec2(-0.4,-0.),0.35),0.0);

  return unionR(s1,unionR(subtractOp(s2,s3),intersectOp(s4,s5)));
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

vec2 DF_vector(vec2 p)
{
	float x = (sceneSDF(vec2(p.x + EPSILON, p.y)).x - sceneSDF(vec2(p.x - EPSILON, p.y)).x)*(0.5 / EPSILON);
	float y = (sceneSDF(vec2(p.x, p.y + EPSILON)).x - sceneSDF(vec2(p.x, p.y - EPSILON)).x)*(0.5 / EPSILON);
	return vec2(x, y);
}

vec3 ray_march(in vec2 org,in vec2 odir) {
    float t = 0.0;
    vec2 dir = odir;
    vec3 total = vec3(0.0);
    vec2 offset = vec2(0.0);
    vec2 r = sceneSDF(org+offset);
    for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
        float dist = abs(r.x);
        if (dist < EPSILON) {
             vec2 n = DF_vector(org+offset);
             float cosq = dot(n,dir);

             float refrac = cosq > 0.0 ? 1.0/1.2:1.2;
             n = cosq > 0.0 ? n : -n;
             if(!is_total_reflection(dir, n, refrac))
             {
                dir = refractR(dir, n, refrac);
                offset += dir*EPSILON*10.0;
                total += vec3(r.y);
             }
			       else return total+vec3(r.y);
        }
        offset += dist*dir;
        t += dist;
        if (t >= MAX_DIST) {
            return total;
        }
        r = sceneSDF(org+offset);
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
        total += ray_march(p, dir);
        seed = mod( seed*1.1234567893490423, 13. );
    }
    total/=float(SAMPLES);
    total = pow(clamp(total,0.0,1.0),vec3(0.45));
    fragColor = vec4(total, 1.0);
}
