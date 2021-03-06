
#version 430

#ifdef COMPUTE_SHADER

#define M_PI 3.1415926535897932384626433832795

struct Triangle
{
    vec3 a;		// sommet
    //int pada;
    vec3 ab;	// arete 1
    //int padb;
    vec3 ac;	// arete 2
    int id;
};

struct BBox
{
    vec3 pmin;
    vec3 pmax;
};

struct NodeGPU
{
    vec3 pmin;
    int left;
    vec3 pmax;
    int right;
};

struct NodeGPU_cousu
{
    vec3 pmin;
    int left;
    vec3 pmax;
    int right;
    vec3 pad;
    int skip;
};

bool internal(in NodeGPU node ) { return (node.right > 0); };                        // renvoie vrai si le noeud est un noeud interne
int internal_left(in NodeGPU node ) {  return node.left; };     // renvoie le fils gauche du noeud interne 
int internal_right(in NodeGPU node ) { return node.right; };   // renvoie le fils droit

bool leaf(in NodeGPU node ) { return (node.right < 0); };                            // renvoie vrai si le noeud est une feuille
int leaf_begin(in NodeGPU node ) { return -node.left; };           // renvoie le premier objet de la feuille
int leaf_end(in NodeGPU node ) { return -node.right; };   


bool leaf(in NodeGPU_cousu node ) { return (node.right < 0); };                            // renvoie vrai si le noeud est une feuille
int leaf_begin(in NodeGPU_cousu node ) { return -node.left; };           // renvoie le premier objet de la feuille
int leaf_end(in NodeGPU_cousu node ) { return -node.right; };   
// int leaf_begin(in NodeGPU_cousu node ) { return (-node.left >> 4 ); };           // renvoie le premier objet de la feuille
// int leaf_end(in NodeGPU_cousu node ) { return ((-node.left >> 4) + (-node.left & 0xF) ); };   


vec3 global( const in vec3 n) { 
    
    float sign= n.z < 0 ? -1.0f : 1.0f;
    float a= -1.0f / (sign + n.z);
    float d= n.x * n.y * a;
    vec3 t= vec3(1.0f + sign * n.x * n.x * a, sign * d, -sign * n.x);
    vec3 b= vec3(d, sign + n.y * n.y * a, -n.y);
    return  vec3(n.x * t + n.y * b + n.z * n); 
};


const uint rng_a= 1103515245;
const uint rng_b= 12345;
const uint rng_mask= (1u << 31) -1u;
const float rng_m= 1u << 31;

// renvoie un reel aleatoire dans [0 1]
float rng( inout uint state )
{
    state= (rng_a * state + rng_b) % uint(rng_m);
    return float(state) / rng_m;
};


// genere une direction sur l'hemisphere,
// cf GI compendium, eq 35
vec3 sample35(in const float u1, in const float u2)
{
    // coordonnees theta, phi
    float cos_theta = sqrt(u1);
    float phi = float(2 * M_PI) * u2;

    // passage vers x, y, z
    float sin_theta = sqrt(1 - cos_theta * cos_theta);
    return vec3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
};

// evalue la densite de proba, la pdf de la direction, cf GI compendium, eq 35
float pdf35(in const vec3 w)
{
    if (w.z < 0)
        return 0;
    return w.z / float(M_PI);
};


// shader storage buffer 0
layout(std430, binding= 0) readonly buffer triangleData
{
    Triangle triangles[];
};

layout(std430, binding= 1) readonly buffer nodeData
{
    NodeGPU nodes[];
};

layout(std430, binding= 2) readonly buffer nodeData_cousu
{
    NodeGPU_cousu nodes_cousu[];
};

BBox bounds( in const NodeGPU node ) 
{
   BBox box;
   box.pmin= node.pmin;
   box.pmax= node.pmax;
   return box;
};

BBox bounds( in const NodeGPU_cousu node ) 
{
   BBox box;
   box.pmin= node.pmin;
   box.pmax= node.pmax;
   return box;
};


struct RayHit
{
    vec3 o;            // origine
    float t;            // p(t)= o + td, position du point d'intersection sur le rayon
    vec3 d;           // direction
    int triangle_id;    // indice du triangle dans le mesh
    vec3 tr_ab;  //cote du trianlge intersecté pour obtenir la normale
    vec3 tr_ac;
    float u, v;

};



bool intersect_new( in const Triangle triangle, inout RayHit ray, in const float tmax)
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
    ray.tr_ab = triangle.ab;
    ray.tr_ac = triangle.ac;

    return true;//(ray.t < tmax && ray.t > 0);
};


bool intersect( const in BBox node, in const RayHit ray, const in vec3 invd)
{
    vec3 dmin= (node.pmin - ray.o) * invd;
    vec3 dmax= (node.pmax - ray.o) * invd;
    
    vec3 tmin= min(dmin, dmax);
    vec3 tmax= max(dmin, dmax);
    float t0= max(tmin.x, max(tmin.y, tmin.z));
    float t1= min(tmax.x, min(tmax.y, tmax.z));

    return max(t0, 0) <= min(t1, 1000.0);
};

uniform mat4 invMatrix;
uniform int frame;
uniform int root;



bool intersect( in BBox bbox, in const RayHit ray )
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    return intersect(bbox, ray, invd);
};

void intersect(inout RayHit ray, in const vec3 invd, in const float tmax)
    {

    int stack[16];
    int top= 0;
    
    // empiler la racine
    stack[top++]= root;//root;
    
    //float tmax= ray.tmax;
    // tant qu'il y a un noeud dans la pile
    while(top > 0)
    {
        int index= stack[--top];
        
        //const NodeGPU node= nodes[index];
        const NodeGPU node= nodes[index];
        //if(intersect(node.bounds, ray, invd))

        if(leaf(node))
        {
            for(int i= leaf_begin(node); i < leaf_end(node); i++)
            {
                intersect_new(triangles[i],ray,tmax);
            }
        }
        else // if(node.internal())
        {
            if(intersect(bounds(node), ray, invd))
            {
                // le noeud est touche, empiler les fils
                stack[top++]= internal_left(node);
                stack[top++]= internal_right(node);
            }
        }
        
    }
};

// parcours iterarif, arbre cousu
// next() : renvoie le prochain noeud dans le sous arbre
// skip( ) : renvoie le prochain sous arbre
void intersect_cousu( inout RayHit ray, in const vec3 invd, in const float tmax)
{

    int index= root;
    while(index != -1)
    {
        const NodeGPU_cousu node= nodes_cousu[index];
        
        if(leaf(node))
        {
            for(int i= leaf_begin(node); i < leaf_end(node); i++){
                intersect_new(triangles[i],ray,tmax);
            }
            index= node.skip;  // passer au prochain sous arbre
        }
            
        else
        {
            if(intersect(bounds(node), ray, invd))
                index= index-1;     // parcourir le sous arbre
            else
                index= node.skip;     // passer au prochain sous arbre
        }
    }
};



void intersect_cousu( inout RayHit ray , in const float tmax) //avec bvh cousu
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    intersect_cousu(ray, invd, tmax);

};


void intersect( inout RayHit ray , in const float tmax) //avec bvh pile
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    intersect(ray, invd, tmax);

};



void test_intersect( inout RayHit rayhit , in const float tmax) //sans bvh
{
    for(int i= 0; i < triangles.length(); i++)
    {
        intersect_new(triangles[i], rayhit, tmax);
    }
};


bool test_intersect_bool( inout RayHit rayhit , in const float tmax) 
{
    for(int i= 0; i < triangles.length(); i++)
    {
        if(intersect_new(triangles[i], rayhit, tmax))
        {
            return true;
        }
    }
    return false;
};



bool intersect_bool(inout RayHit ray, in const vec3 invd, in const float tmax)
{

    int stack[16];
    int top= 0;
    
    // empiler la racine
    stack[top++]= root;//root;
    
    //float tmax= ray.tmax;
    // tant qu'il y a un noeud dans la pile
    while(top > 0)
    {
        int index= stack[--top];
        const NodeGPU node= nodes[index];


        if(leaf(node))
        {                
 
            for(int i= leaf_begin(node); i < leaf_end(node); i++){
                if(intersect_new(triangles[i],ray,tmax)){
                    return true;
                };
            }
        }
        else // if(node.internal())
        {
            if(intersect(bounds(node), ray, invd))
            {
            stack[top++]= internal_left(node);
            stack[top++]= internal_right(node);
            }
        }
        
    }
    return false;
};


bool intersect_bool( inout RayHit ray , in const float tmax)
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    return intersect_bool(ray, invd, tmax);

};

bool intersect_bool_cousu( inout RayHit ray, in const vec3 invd, in const float tmax)
{

    int index= root;
    while(index != -1)
    {
        const NodeGPU_cousu node= nodes_cousu[index];
        
        if(leaf(node))
        {
            for(int i= leaf_begin(node); i < leaf_end(node); i++){
                if(intersect_new(triangles[i],ray,tmax)){
                    return true;
                };
            }
            index= node.skip;  // passer au prochain sous arbre
        }
            
        else
        {
            if(intersect(bounds(node), ray, invd))
                index= index-1;     // parcourir le sous arbre
            else
                index= node.skip;     // passer au prochain sous arbre
        }
    }
    return false;
};


bool intersect_bool_cousu( inout RayHit ray , in const float tmax)
{
    vec3 invd= vec3(1.0f / ray.d.x, 1.0f / ray.d.y, 1.0f / ray.d.z);
    return intersect_bool_cousu(ray, invd, tmax);

};

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

    float tmax =1000.0;
    float hit= tmax;	// tmax = far, une intersection valide est plus proche que l'extremite du rayon / far...
    float hitu= 0.0;
    float hitv= 0.0;
    int id;
        //float t, u, v;    
    float t, u, v;  


    RayHit rayhit = RayHit(o,hit,d,id,vec3(-1),vec3(-1),u,v);
    

    //////////////////////// sans bvh :  /////////////////////////////////
    //test_intersect( rayhit , tmax) ;

    //////////////////////// avec bvh pile:  /////////////////////////////////
    //intersect( rayhit , tmax) ;

    //////////////////////// avec bvh cousu:  /////////////////////////////////
    intersect_cousu( rayhit , tmax) ;

    hit=rayhit.t;
    hitu=rayhit.u;
    hitv=rayhit.v;
    id=rayhit.triangle_id;
    
    float w = 1.0 - hitu - hitv;
    vec3 p = rayhit.o+ rayhit.t * rayhit.d;
    vec3 n_p = normalize(cross(rayhit.tr_ab,rayhit.tr_ac));

    int N_ray = 4;
    vec4 ambient = vec4(0.0);
    //uint state = 0;

    uint state = imageLoad(seeds, ivec2(gl_GlobalInvocationID.xy)).x;

    
    for (int ni = 0; ni<N_ray; ni++){

        uint k =  uint((frame*N_ray+ni));

        float u1 = rng(state);
        float u2 = rng(state);
        
        vec3 d_l = sample35(u1,u2);
        float pdf = pdf35(d_l);


        //passage dans le repere monde
        float sign= n_p.z < 0 ? -1.0f : 1.0f;
        float a= -1.0f / (sign + n_p.z);
        float e= n_p.x * n_p.y * a;
        vec3 t= vec3(1.0f + sign * n_p.x * n_p.x * a, sign * e, -sign * n_p.x);
        vec3 b= vec3(e, sign + n_p.y * n_p.y * a, -n_p.y);
        d_l = d_l.x * t + d_l.y * b + d_l.z * n_p; 
        
        
        //d_l = global(d_l);
        int v = 0;
        float vu_sun = 1;
        

        float tmax2=1000.0;
        float hit2= tmax2;	// tmax = far, une intersection valide est plus proche que l'extremite du rayon / far...
        float hitu2= 0.0;
        float hitv2= 0.0;
        int id2=-1;
        float t2, u22, v2;
        //int i2;
        RayHit rayhit2 = RayHit(p+0.01*n_p,hit2,d_l,id2,vec3(-1),vec3(-1),u22,v2);


        //////////////////////// sans bvh :  /////////////////////////////////
        // if(test_intersect_bool(rayhit2 , tmax2)){
        //     vu_sun = 0;
        // }

        //////////////////////// avec bvh pile:  /////////////////////////////////
        // if(intersect_bool(rayhit2 , tmax2)){
        //     vu_sun = 0;
        // }

        //////////////////////// avec bvh cousu:  /////////////////////////////////
        if(intersect_bool_cousu(rayhit2 , tmax2)){
            vu_sun = 0;
        }


        if (vu_sun==1){

            vec4 diffuse = vec4(1.0);
            float cos_theta = max(float(0.0), dot(n_p, d_l));

            ambient = ambient + diffuse*cos_theta / (float(M_PI)*pdf);
        }

    }

    vec4 Color = vec4(0.0);
    if (hit < 1.0 && hit > 0){ 
        Color = Color + ambient;
    }
    vec4 curr_im = imageLoad(image, ivec2(gl_GlobalInvocationID.xy));

    Color = (curr_im *frame*N_ray + Color) / (frame *N_ray+N_ray);

    // ecrire le resultat dans l'image
    imageStore(image, ivec2(gl_GlobalInvocationID.xy), Color );
    imageStore(seeds, ivec2(gl_GlobalInvocationID.xy), uvec4(state) );
};
#endif
