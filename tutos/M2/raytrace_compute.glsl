
#version 430

#ifdef COMPUTE_SHADER

struct Triangle
{
    vec3 a;		// sommet
    vec3 ab;	// arete 1
    vec3 ac;	// arete 2
    int id;
};

struct BBox
{
    vec3 pmin;
    int paspmin;
    vec3 pmax;
    int paspmax;
};


struct Node
{
    BBox bounds;
    int left;
    int right;
    int padleft;
    int padright;
};

bool internal( Node node ) { return (node.right > 0); };                        // renvoie vrai si le noeud est un noeud interne
int internal_left( Node node ) {  return node.left; };     // renvoie le fils gauche du noeud interne 
int internal_right( Node node ) { return node.right; };   // renvoie le fils droit

bool leaf( Node node ) { return (node.right < 0); };                            // renvoie vrai si le noeud est une feuille
int leaf_begin( Node node ) { return -node.left; };           // renvoie le premier objet de la feuille
int leaf_end( Node node ) { return -node.right; };   


vec3 global( const in vec3 n) { 
    
    float sign= n.z < 0 ? -1.0f : 1.0f;
    float a= -1.0f / (sign + n.z);
    float d= n.x * n.y * a;
    vec3 t= vec3(1.0f + sign * n.x * n.x * a, sign * d, -sign * n.x);
    vec3 b= vec3(d, sign + n.y * n.y * a, -n.y);
    return  vec3(n.x * t + n.y * b + n.z * n); 
};

/*
struct World
{
    World( const vec3 _n ) : n(_n){
        float sign= std::copysign(1.0f, n.z);
        float a= -1.0f / (sign + n.z);
        float d= n.x * n.y * a;
        t= Vector(1.0f + sign * n.x * n.x * a, sign * d, -sign * n.x);
        b= Vector(d, sign + n.y * n.y * a, -n.y);
    }

    // transforme le vecteur du repere local vers le repere du monde
    Vector operator( ) ( const Vector& local )  const { return local.x * t + local.y * b + local.z * n; }

    // transforme le vecteur du repere du monde vers le repere local
    Vector inverse( const Vector& global ) const { return Vector(dot(global, t), dot(global, b), dot(global, n)); }

    Vector t;
    Vector b;
    Vector n;
};

*/

const uint rng_a= 1103515245;
const uint rng_b= 12345;
const uint rng_mask= (1u << 31) -1u;
const float rng_m= 1u << 31;

// renvoie un reel aleatoire dans [0 1]
float rng( inout uint state )
{
    state= (rng_a * state + rng_b) % uint(rng_m);
    return float(state) / rng_m;
}

#define M_PI 3.1415926535897932384626433832795


// genere une direction sur l'hemisphere,
// cf GI compendium, eq 35
vec3 sample35(const float u1, const float u2)
{
    // coordonnees theta, phi
    float cos_theta = sqrt(u1);
    float phi = float(2 * M_PI) * u2;

    // passage vers x, y, z
    float sin_theta = sqrt(1 - cos_theta * cos_theta);
    return vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
}

// evalue la densite de proba, la pdf de la direction, cf GI compendium, eq 35
float pdf35(const vec3 w)
{
    if (w.z < 0)
        return 0;
    return w.z / float(M_PI);
}


// shader storage buffer 0
layout(std430, binding= 0) readonly buffer triangleData
{
    Triangle triangles[];
};

layout(std430, binding= 1) readonly buffer nodeData
{
    Node nodes[];
};


// bool intersect( const Triangle triangle, const vec3 o, const vec3 d, const float tmax, out float rt,/* out int id, */out float ru, out float rv )
// {
//     vec3 pvec= cross(d, triangle.ac);
//     float det= dot(triangle.ab, pvec);
//     float inv_det= 1.0f / det;
    
//     vec3 tvec= o - triangle.a;
//     float u= dot(tvec, pvec) * inv_det;
//     vec3 qvec= cross(tvec, triangle.ab);
//     float v= dot(d, qvec) * inv_det;
    
//     /* calculate t, ray intersects triangle */
//     rt= dot(triangle.ac, qvec) * inv_det;
//     ru= u;
//     rv= v;
//     //id = triangle.id;
    
//     // ne renvoie vrai que si l'intersection est valide : 
//     // interieur du triangle, 0 < u < 1, 0 < v < 1, 0 < u+v < 1
//     if(any(greaterThan(vec3(u, v, u+v), vec3(1, 1, 1))) || any(lessThan(vec2(u, v), vec2(0, 0))))
//         return false;
//     // comprise entre 0 et tmax du rayon
//     return (rt < tmax && rt > 0);
// };

struct RayHit
{
    vec3 o;            // origine
    float t;            // p(t)= o + td, position du point d'intersection sur le rayon
    vec3 d;           // direction
    int triangle_id;    // indice du triangle dans le mesh
    float u, v;

};

// bool intersect( const Triangle triangle, const vec3 o, const vec3 d, const float tmax, out float rt, out int id, out float ru, out float rv, RayHit ray )
// {
//     vec3 pvec= cross(ray.d, triangle.ac);
//     float det= dot(triangle.ab, pvec);
//     float inv_det= 1.0f / det;
    
//     vec3 tvec= ray.o - triangle.a;
//     float u= dot(tvec, pvec) * inv_det;
//     vec3 qvec= cross(tvec, triangle.ab);
//     float v= dot(d, qvec) * inv_det;
    
//     ray.t= dot(triangle.ac, qvec) * inv_det;
//     ray.triangle_id = triangle.id;
//     ray.u= u;
//     ray.v= v;

//     /* calculate t, ray intersects triangle */
//     // rt= dot(triangle.ac, qvec) * inv_det;
//     // ru= u;
//     // rv= v;
//     // id = triangle.id;
    
//     // ne renvoie vrai que si l'intersection est valide : 
//     // interieur du triangle, 0 < u < 1, 0 < v < 1, 0 < u+v < 1
//     if(any(greaterThan(vec3(u, v, u+v), vec3(1, 1, 1))) || any(lessThan(vec2(u, v), vec2(0, 0))))
//         return false;
//     // comprise entre 0 et tmax du rayon
//     return (rt < tmax && rt > 0);
// };


bool intersect_new( const Triangle triangle, inout RayHit ray, const float tmax)
{
    
    //ray = ray;
    vec3 pvec= cross(ray.d, triangle.ac);
    float det= dot(triangle.ab, pvec);
    float inv_det= 1.0f / det;
    
    vec3 tvec= ray.o - triangle.a;
    float u= dot(tvec, pvec) * inv_det;
    if(u < 0.0 || u > 1.0) return false;
    
    vec3 qvec= cross(tvec, triangle.ab);
    float v= dot(ray.d, qvec) * inv_det;
    if(v < 0.0 || u + v > 1.0) return false;
    
    float t= dot(triangle.ac, qvec) * inv_det;
    if(t < 0.0 || t > ray.t) return false;



    ray.t = t;
    ray.triangle_id = triangle.id;
    ray.u= u;
    ray.v= v;

    /* calculate t, ray intersects triangle */
    // rt= dot(triangle.ac, qvec) * inv_det;
    // ru= u;
    // rv= v;
    // id = triangle.id;


    //return true;
    
    // ne renvoie vrai que si l'intersection est valide : 
    // interieur du triangle, 0 < u < 1, 0 < v < 1, 0 < u+v < 1
    if(any(greaterThan(vec3(u, v, u+v), vec3(1, 1, 1))) || any(lessThan(vec2(u, v), vec2(0, 0)))){
        return false;
    }
    // comprise entre 0 et tmax du rayon
    return (ray.t < tmax && ray.t > 0);
};

bool intersect( const Triangle triangle, inout RayHit ray, const float tmax)
{
    
    //ray = ray;
    vec3 pvec= cross(ray.d, triangle.ac);
    float det= dot(triangle.ab, pvec);
    float inv_det= 1.0f / det;
    
    vec3 tvec= ray.o - triangle.a;
    float u= dot(tvec, pvec) * inv_det;
    if(u < 0.0 || u > 1.0) return false;
    
    vec3 qvec= cross(tvec, triangle.ab);
    float v= dot(ray.d, qvec) * inv_det;
    if(v < 0.0 || u + v > 1.0) return false;
    
    float t= dot(triangle.ac, qvec) * inv_det;
    if(t < 0.0 || t > ray.t) return false;



    ray.t = t;
    ray.triangle_id = triangle.id;
    ray.u= u;
    ray.v= v;

    /* calculate t, ray intersects triangle */
    // rt= dot(triangle.ac, qvec) * inv_det;
    // ru= u;
    // rv= v;
    // id = triangle.id;


    //return true;
    
    // ne renvoie vrai que si l'intersection est valide : 
    // interieur du triangle, 0 < u < 1, 0 < v < 1, 0 < u+v < 1
    if(any(greaterThan(vec3(u, v, u+v), vec3(1, 1, 1))) || any(lessThan(vec2(u, v), vec2(0, 0)))){
        return false;
    }
    // comprise entre 0 et tmax du rayon
    return (ray.t < tmax && ray.t > 0);
};


bool intersect2( const Triangle triangle, inout RayHit ray, const float tmax)
{
    vec3 pvec= cross(ray.d, triangle.ac);
    float det= dot(triangle.ab, pvec);
    float inv_det= 1.0f / det;
    
    vec3 tvec= ray.o - triangle.a;
    float u= dot(tvec, pvec) * inv_det;
    vec3 qvec= cross(tvec, triangle.ab);
    float v= dot(ray.d, qvec) * inv_det;
    
    ray.t= dot(triangle.ac, qvec) * inv_det;
    ray.triangle_id = triangle.id;
    ray.u= u;
    ray.v= v;

    /* calculate t, ray intersects triangle */
    // rt= dot(triangle.ac, qvec) * inv_det;
    // ru= u;
    // rv= v;
    // id = triangle.id;


    //return true;
    
    // ne renvoie vrai que si l'intersection est valide : 
    // interieur du triangle, 0 < u < 1, 0 < v < 1, 0 < u+v < 1
    if(any(greaterThan(vec3(u, v, u+v), vec3(1, 1, 1))) || any(lessThan(vec2(u, v), vec2(0, 0)))){
        return false;
    }
    // comprise entre 0 et tmax du rayon
    return (ray.t < tmax && ray.t > 0);
};




void swap(inout float a, inout float b){
    float tempa = a;
    a = b;
    b = tempa;
}

bool intersect( BBox bbox, const RayHit ray, const vec3 invd )
{
    vec3 rmin= bbox.pmin;
    vec3 rmax= bbox.pmax;
    if(ray.d.x < 0) {
        float temp = rmin.x;
        rmin.x = rmax.x;
        rmax.x = temp;
        //std::swap(rmin.x, rmax.x);
    }
    if(ray.d.y < 0){
        float temp = rmin.y;
        rmin.y = rmax.y;
        rmax.y = temp;    
        //std::swap(rmin.y, rmax.y);
    } 
    if(ray.d.z < 0){
        float temp = rmin.z;
        rmin.z = rmax.z;
        rmax.z = temp;   
        //std::swap(rmin.z, rmax.z);
    } 
    vec3 dmin= (rmin - ray.o) * invd;
    vec3 dmax= (rmax - ray.o) * invd;
    
    float tmin= max(dmin.z, max(dmin.y, max(dmin.x, 0.f)));
    float tmax= min(dmax.z, min(dmax.y, min(dmax.x, ray.t)));
    return tmin<= tmax;
};


bool intersect( BBox bbox, const RayHit ray )
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    return intersect(bbox, ray, invd);
};

void intersect(inout RayHit ray, const vec3 invd, const float tmax)
    {

    int stack[64];
    int top= 0;
    
    // empiler la racine
    stack[top++]= 64;//root;
    
    //float tmax= ray.tmax;
    // tant qu'il y a un noeud dans la pile
    while(top > 0)
    {
        int index= stack[--top];
        
        const Node node= nodes[index];
        if(intersect(node.bounds, ray, invd))
        {
            if(leaf(node))
            {
                for(int i= leaf_begin(node); i < leaf_end(node); i++){
                    intersect_new(triangles[i],ray,tmax);
                }
            }
            else // if(node.internal())
            {
                //assert(top +1 < 64);       // le noeud est touche, empiler les fils
                stack[top++]= internal_left(node);
                stack[top++]= internal_right(node);
            }
        }
    }
};


void intersect( RayHit ray , const float tmax)
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    intersect(ray, invd, tmax);

};


void test_intersect( inout RayHit rayhit , const float tmax) 
{
    for(int i= 0; i < triangles.length(); i++)
    {
        intersect_new(triangles[i], rayhit, tmax);
    }
};

void test_intersect_old( inout RayHit rayhit , const float tmax, inout float hit, inout float hitu, inout float hitv, inout int id) 
{
    for(int i= 0; i < triangles.length(); i++)
    {
        //float t, u, v;    
        // float t, u, v;  
        // //int i2;
        // int i2;
        // RayHit rayhit = RayHit(o,t,d,i2,u,v);
        
        //if(intersect(triangles[i], o, d, hit, t, u, v))
        if(intersect2(triangles[i], rayhit, tmax ))
        //out float rt, out int id, out float ru, out float rv)
        //if(intersect(triangles[i], o, d, hit, t, i2, u, v, rayhit))
        {
            if (hit>rayhit.t){ //zbuffer
                hit= rayhit.t;
                hitu=rayhit.u;//1-t;//
                hitv=rayhit.v; //1-t;//
                id=rayhit.triangle_id;//rayhit.triangle_id;//triangles[i].id;
            }
        }
    }
};

bool test_intersect_bool( inout RayHit rayhit , const float tmax, inout float hit, inout float hitu, inout float hitv, inout int id) 
{
    for(int i= 0; i < triangles.length(); i++)
    {
        //float t, u, v;    
        // float t, u, v;  
        // //int i2;
        // int i2;
        // RayHit rayhit = RayHit(o,t,d,i2,u,v);
        
        //if(intersect(triangles[i], o, d, hit, t, u, v))
        if(intersect2(triangles[i], rayhit, tmax))
        //out float rt, out int id, out float ru, out float rv)
        //if(intersect(triangles[i], o, d, hit, t, i2, u, v, rayhit))
        {
            hit= rayhit.t;
            hitu=rayhit.u;//1-t;//
            hitv=rayhit.v; //1-t;//
            id=rayhit.triangle_id;//rayhit.triangle_id;//triangles[i].id;
            return true;
        }
    }
    return false;
};


uniform mat4 invMatrix;
uniform int frame;

// image resultat
layout(binding= 0, rgba8)  coherent uniform image2D image;
layout(binding= 1, r32ui)  coherent uniform uimage2D seeds;


// 8x8 threads
layout( local_size_x= 8, local_size_y= 8 ) in;
void main( )
{
    // recupere le threadID 2d, et l'utilise directement comme coordonnees de pixel
    vec2 position= vec2(gl_GlobalInvocationID.xy);
    
    // construction du rayon pour le pixel, passage depuis le repere image vers le repere monde
    vec4 oh= invMatrix * vec4(position, 0, 1);       // origine sur near
    vec4 eh= invMatrix * vec4(position, 1, 1);       // extremite sur far

    // origine et direction
    vec3 o= oh.xyz / oh.w;                              // origine
    vec3 d= eh.xyz / eh.w - oh.xyz / oh.w;              // direction

    float tmax =1.0;
    float hit= 1.0;	// tmax = far, une intersection valide est plus proche que l'extremite du rayon / far...
    float hitu= 0.0;
    float hitv= 0.0;
    int id;
        //float t, u, v;    
    float t, u, v;  
    //int i2;
    int i2;
    RayHit rayhit = RayHit(o,hit,d,i2,u,v);
    //intersect( rayhit , tmax, hit, hitu, hitv,id);
    test_intersect( rayhit , tmax) ;

    // if (rayhit.t < 1.0 && rayhit.t > 0){ //if dans mesh
    hit= rayhit.t;
    hitu=rayhit.u;//1-t;//
    hitv=rayhit.v; //1-t;//
    id=rayhit.triangle_id;//rayhit.triangle_id;//triangles[i].id;
    //}


    // for(int i= 0; i < triangles.length(); i++)
    // {

    //     //if(intersect(triangles[i], o, d, hit, t, u, v))
    //     if(intersect(triangles[i], rayhit, hit))
    //     //out float rt, out int id, out float ru, out float rv)
    //     //if(intersect(triangles[i], o, d, hit, t, i2, u, v, rayhit))
    //     {
    //         hit=rayhit.t;
    //         hitu=rayhit.u;//1-t;//
    //         hitv=rayhit.v; //1-t;//
    //         id=rayhit.triangle_id;//rayhit.triangle_id;//triangles[i].id;
    //     }
    // }
    
    float w = 1.0 - hitu - hitv;
    vec3 p = triangles[id].a + hitu*triangles[id].ab + hitv*triangles[id].ac;
    vec3 n_p = normalize(cross(triangles[id].ab,triangles[id].ac));

    int N_ray = 8;
    vec4 ambient = vec4(0.0);
    //uint state = 0;

    uint state = imageLoad(seeds, ivec2(gl_GlobalInvocationID.xy)).x;

    
    for (int ni = 0; ni<N_ray; ni++){


 
        uint k = uint(frame*N_ray+state);

        //vec3 d_l = normalize(vec3(0.0,0.0,1.98)-p);///
        //vec3 d_l = vec3(0.0,1.0,0.0);
        float u1 = rng(k);
        float u2 = rng(k);
        
        vec3 d_l = sample35(u1,u2);
        float pdf = pdf35(d_l);



        // if(dot(d_l, n_p) > 0){
        //     n_p= -n_p;
        // }

        float sign= n_p.z < 0 ? -1.0f : 1.0f;
        float a= -1.0f / (sign + n_p.z);
        float e= n_p.x * n_p.y * a;
        vec3 t= vec3(1.0f + sign * n_p.x * n_p.x * a, sign * e, -sign * n_p.x);
        vec3 b= vec3(e, sign + n_p.y * n_p.y * a, -n_p.y);
        d_l = d_l.x * t + d_l.y * b + d_l.z * n_p; 
        
        
        //d_l = global(d_l);
        int v = 0;
        float vu_sun = 1;
        


        float hit2= 1.0;	// tmax = far, une intersection valide est plus proche que l'extremite du rayon / far...
        float hitu2= 0.0;
        float hitv2= 0.0;
        int id2;
        float t2, u22, v2;
        //int i2;
        RayHit rayhit2 = RayHit(p+0.01*n_p,t2,d_l,id2,u22,v2);

        if(test_intersect_bool(rayhit2 , 10000, hit2, hitu2, hitv2, id2)==true){
            v = 0;
            vu_sun = 0;
        }

        // for(int j= 0; j < triangles.length(); j++)
        // {

        //     //if(intersect(triangles[j], p+0.01*n_p, d_l, 1000, t2, /*id2,**/ u2, v2))
        //     if(intersect(triangles[j], rayhit2, 1000))
        //     {
        //         // hit2= t2;
        //         // hitu2=u2;//1-t;//
        //         // hitv2=v2; //1-t;//
        //         // id2=j;//rayhit.triangle_id;
        //         v = 0;
        //         vu_sun = 0;
        //         break;
                
        //     }
        // }

        if (vu_sun==1){
            //vec4 diffuse = vec4(0,0.140,0.450,0.0911);
            vec4 diffuse = vec4(1.0);
            float cos_theta = max(float(0.0), dot(n_p, d_l));
            //float cos_theta = abs(dot(n_p, d_l));
            ambient = ambient + diffuse*cos_theta / (float(M_PI)*pdf);// * (cos_theta / (float(M_PI) * pdf35 * N_ray));
        }

    }
    // vec4 Color = vec4(0, 0, 0, 1);

    // if (hit < 1.0 && hit > 0){ //if dans mesh
    //     //Color = vec4(hitu, hitv, 0, 1)*v;
    //     Color += vec4(1-hit, 1-hit, 1-hit, 1)*(1-v);
    //     Color += vec4(1-hit, 1-hit, 0, 1)*v;
    // }
    // else {
    //     Color = vec4(0, 0, 1, 1);
    // }
    vec4 Color = vec4(0.0);
    if (hit < 1.0 && hit > 0){ //if dans mesh
        Color = Color + ambient;//vec4(0, 0, 0, 1);
    }
    vec4 curr_im = imageLoad(image, ivec2(gl_GlobalInvocationID.xy));
    Color = (curr_im *frame*N_ray + Color) / (frame *N_ray+N_ray);
    // ecrire le resultat dans l'image
    imageStore(image, ivec2(gl_GlobalInvocationID.xy), Color );
}
#endif
