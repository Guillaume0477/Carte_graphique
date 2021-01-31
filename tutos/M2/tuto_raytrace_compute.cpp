
#include <cfloat>
#include <cmath>
#include <algorithm>

#include "app_time.h"

#include "vec.h"
#include "color.h"
#include "mat.h"

#include "mesh.h"
#include "wavefront.h"

#include "program.h"
#include "uniforms.h"

#include "orbiter.h"
#include <random>


// cf tuto_storage

namespace glsl 
{
    // type de base alignes sur 4 octets
    template < typename T >
    struct alignas(4) gscalar
    { 
        alignas(4) T x;
        
        gscalar( ) : x(T()) {}
        gscalar( const T& v ) : x(v) {}
        gscalar& operator= ( const T& v ) { x= v; return *this; }
        operator T ( ) { return x; }
    };

    typedef gscalar<float> gfloat;
    typedef gscalar<int> gint;
    typedef gscalar<unsigned int> guint;
    typedef gscalar<bool> gbool;
    
    // vec2, alignes sur 2 * alignement type de base du vecteur
    template < typename T >
    struct alignas(8) gvec2
    {
        alignas(4) T x, y;
        
        gvec2( ) {}
        gvec2( const gvec2<T>& v ) : x(v.x), y(v.y) {}
        gvec2( const ::vec2& v ) : x(v.x), y(v.y) {}
        gvec2& operator= ( const gvec2<T>& v ) { x= v.x; y= v.y; return *this; }
        gvec2& operator= ( const ::vec2& v ) { x= v.x; y= v.y; return *this; }
        operator ::vec2 ( ) { return ::vec2(float(x), float(y)); }
    };
    
    typedef gvec2<float> vec2;
    typedef gvec2<int> ivec2;
    typedef gvec2<unsigned int> uvec2;
    typedef gvec2<int> bvec2;
    
    // vec3, alignes sur 4 * alignement type de base du vecteur
    template < typename T >
    struct alignas(16) gvec3
    {
        alignas(4) T x, y, z;
        
        gvec3( ) {}
        gvec3( const gvec3<T>& v ) : x(v.x), y(v.y), z(v.z) {}
        gvec3( const ::vec3& v ) : x(v.x), y(v.y), z(v.z) {}
        gvec3( const Point& v ) : x(v.x), y(v.y), z(v.z) {}
        gvec3( const Vector& v ) : x(v.x), y(v.y), z(v.z) {}
        gvec3& operator= ( const gvec3<T>& v ) { x= v.x; y= v.y; z= v.z; return *this; }
        gvec3& operator= ( const ::vec3& v ) { x= v.x; y= v.y; z= v.z; return *this; }
        gvec3& operator= ( const Point& v ) { x= v.x; y= v.y; z= v.z; return *this; }
        gvec3& operator= ( const Vector& v ) { x= v.x; y= v.y; z= v.z; return *this; }
        operator ::vec3 ( ) { return ::vec3(float(x), float(y), float(y)); }
    };
    
    typedef gvec3<float> vec3;
    typedef gvec3<int> ivec3;
    typedef gvec3<unsigned int> uvec3;
    typedef gvec3<int> bvec3;
    
    // vec4, alignes sur 4 * alignement type de base du vecteur
    template < typename T >
    struct alignas(16) gvec4
    {
        alignas(4) T x, y, z, w;
        
        gvec4( ) {}
        gvec4( const gvec4<T>& v ) : x(v.x), y(v.y), z(v.z), w(v.w) {}
        gvec4( const ::vec4& v ) : x(v.x), y(v.y), z(v.z), w(v.w) {}
        gvec4& operator= ( const gvec4<T>& v ) { x(v.x), y(v.y), z(v.z), w(v.w) ; return *this; }
        gvec4& operator= ( const ::vec4& v ) { x(v.x), y(v.y), z(v.z), w(v.w) ; return *this; }
        gvec4& operator= ( const Color& c ) { x= c.r; y= c.g; z= c.b; w= c.a; return *this; }
        operator ::vec4 ( ) { return ::vec4(float(x), float(y), float(y), float(w)); }
    };
    
    typedef gvec4<float> vec4;
    typedef gvec4<int> ivec4;
    typedef gvec4<unsigned int> uvec4;
    typedef gvec4<int> bvec4;
}

// struct Triangle 
// {
//     vec3 a;
//     int pada;
//     vec3 ab;
//     int padb;
//     vec3 ac;
//     int id;
//     // int pad1;
//     // int pad2;
//     // int pad3;
// };



struct Ray
{
    Point o;
    float pad;
    Vector d;
    float tmax;

    Ray() : o(), d(), tmax(0) {}
    Ray(const Point &_o, const Point &_e) : o(_o), d(Vector(_o, _e)), tmax(1) {}
    Ray(const Point &_o, const Vector &_d) : o(_o), d(_d), tmax(FLT_MAX) {}
};

// intersection rayon / triangle.
struct Hit
{
    int triangle_id;
    float t;
    float u, v;

    Hit() : triangle_id(-1), t(0), u(0), v(0) {} // pas d'intersection
    Hit(const int _id, const float _t, const float _u, const float _v) : triangle_id(_id), t(_t), u(_u), v(_v) {}

    operator bool() const { return (triangle_id != -1); } // renvoie vrai si l'intersection est initialisee...
};


 struct RayHit
 {
     Point o;            // origine
     float t;            // p(t)= o + td, position du point d'intersection sur le rayon
     Vector d;           // direction
     int triangle_id;    // indice du triangle dans le mesh
     float u, v;
     int x, y;
     
     RayHit( const Point& _o, const Point& _e ) :  o(_o), t(1), d(Vector(_o, _e)), triangle_id(-1), u(), v(), x(), y() {}
     RayHit( const Point& _o, const Point& _e, const int _x, const int _y ) :  o(_o), t(1), d(Vector(_o, _e)), triangle_id(-1), u(), v(), x(_x), y(_y) {}
     operator bool ( ) { return (triangle_id != -1); }
 };
  
  
 struct BBoxHit
 {
     float tmin, tmax;
     
     BBoxHit() : tmin(FLT_MAX), tmax(-FLT_MAX) {}
     BBoxHit( const float _tmin, const float _tmax ) : tmin(_tmin), tmax(_tmax) {}
     float centroid( ) const { return (tmin + tmax) / 2; }
     operator bool( ) const { return tmin <= tmax; }
 };
  


// renvoie la normale interpolee d'un triangle.
Vector normal(const Hit &hit, const TriangleData &triangle)
{
    return normalize((1 - hit.u - hit.v) * Vector(triangle.na) + hit.u * Vector(triangle.nb) + hit.v * Vector(triangle.nc));
}

// renvoie la normale interpolee d'un triangle.
Vector normal(const RayHit &hit, const TriangleData &triangle)
{
    return normalize((1 - hit.u - hit.v) * Vector(triangle.na) + hit.u * Vector(triangle.nb) + hit.v * Vector(triangle.nc));
}

// renvoie le point d'intersection sur le triangle.
Point point(const Hit &hit, const TriangleData &triangle)
{
    return (1 - hit.u - hit.v) * Point(triangle.a) + hit.u * Point(triangle.b) + hit.v * Point(triangle.c);
}

// renvoie le point d'intersection sur le rayon
Point point(const Hit &hit, const Ray &ray)
{
    return ray.o + hit.t * ray.d;
}



struct BBox
 {
     Point pmin;
     int paspmin;
     Point pmax;
     int paspmax;
     
     BBox( ) : pmin(),paspmin(), pmax(),paspmax(){}
     
     BBox( const vec3& p ) : pmin(p),paspmin(1), pmax(p),paspmax(1) {}
     BBox( const BBox& box ) : pmin(box.pmin),paspmin(1), pmax(box.pmax),paspmax(1) {}
     
     BBox& insert( const Point& p ) { pmin= min(pmin, p); pmax= max(pmax, p); return *this; }
     BBox& insert( const BBox& box ) { pmin= min(pmin, box.pmin); pmax= max(pmax, box.pmax); return *this; }
     
     float centroid( const int axis ) const { return (pmin(axis) + pmax(axis)) / 2; }
     
     BBoxHit intersect( const RayHit& ray ) const
     {
         Vector invd= Vector(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
         return intersect(ray, invd);
     }
     
     BBoxHit intersect( const RayHit& ray, const Vector& invd ) const
     {
         Point rmin= pmin;
         Point rmax= pmax;
         if(ray.d.x < 0) std::swap(rmin.x, rmax.x);
         if(ray.d.y < 0) std::swap(rmin.y, rmax.y);
         if(ray.d.z < 0) std::swap(rmin.z, rmax.z);
         Vector dmin= (rmin - ray.o) * invd;
         Vector dmax= (rmax - ray.o) * invd;
         
         float tmin= std::max(dmin.z, std::max(dmin.y, std::max(dmin.x, 0.f)));
         float tmax= std::min(dmax.z, std::min(dmax.y, std::min(dmax.x, ray.t)));
         return BBoxHit(tmin, tmax);
     }
 };




// triangle "intersectable".
struct Triangle
{
    vec3 p;
    int padp;
    vec3 e1;
    int pade1;
    vec3 e2;
    int id;
//     vec3 a;
//     int pada;
//     vec3 ab;
//     int padb;
//     vec3 ac;
//     int id;

    Triangle( const TriangleData& data, const int _id ) : p(data.a),padp(1), e1(Vector(data.a, data.b)),pade1(1), e2(Vector(data.a, data.c)), id(_id) {}
    Triangle( const Point &_a, const int _padp, const Point &_b, const int _pade1 , const Point &_c, const int _id) : p(_a), padp(_padp), e1(Vector(_a, _b)), pade1(_pade1), e2(Vector(_a, _c)), id(_id) {}

    /* calcule l'intersection ray/triangle
        cf "fast, minimum storage ray-triangle intersection" 
        http://www.graphics.cornell.edu/pubs/1997/MT97.pdf
        
        renvoie faux s'il n'y a pas d'intersection valide (une intersection peut exister mais peut ne pas se trouver dans l'intervalle [0 htmax] du rayon.)
        renvoie vrai + les coordonnees barycentriques (u, v) du point d'intersection + sa position le long du rayon (t).
        convention barycentrique : p(u, v)= (1 - u - v) * a + u * b + v * c
    */
    Hit intersect(const Ray &ray, const float htmax) const
    {
        Vector pvec = cross(ray.d, e2);
        float det = dot(e1, pvec);

        float inv_det = 1 / det;
        Vector tvec(p, ray.o);

        float u = dot(tvec, pvec) * inv_det;
        if (u < 0 || u > 1)
            return Hit();

        Vector qvec = cross(tvec, e1);
        float v = dot(ray.d, qvec) * inv_det;
        if (v < 0 || u + v > 1)
            return Hit();

        float t = dot(e2, qvec) * inv_det;
        if (t > htmax || t < 0)
            return Hit();

        return Hit(id, t, u, v); // p(u, v)= (1 - u - v) * a + u * b + v * c
    }
    void intersect( RayHit &ray ) const
    {
        Vector pvec= cross(ray.d, e2);
        float det= dot(e1, pvec);
        
        float inv_det= 1 / det;
        Vector tvec(p, ray.o);
        
        float u= dot(tvec, pvec) * inv_det;
        if(u < 0 || u > 1) return;
        
        Vector qvec= cross(tvec, e1);
        float v= dot(ray.d, qvec) * inv_det;
        if(v < 0 || u + v > 1) return;
        
        float t= dot(e2, qvec) * inv_det;
        if(t < 0 || t > ray.t) return;
        
        // touche !!
        ray.t= t;
        ray.triangle_id= id;
        ray.u= u;
        ray.v= v;
    }
    
    BBox bounds( ) const 
    {
        BBox box(p);
        return box.insert(Point(p)+Vector(e1)).insert(Point(p)+Vector(e2));
    }
};

/*
// ensemble de triangles.
// a remplacer par une vraie structure acceleratrice, un bvh, par exemple
struct BVH
{
    std::vector<Triangle> triangles;

    BVH() = default;
    BVH(const Mesh &mesh) { build(mesh); }

    void build(const Mesh &mesh)
    {
        triangles.clear();
        triangles.reserve(mesh.triangle_count());
        for (int id = 0; id < mesh.triangle_count(); id++)
        {
            TriangleData data = mesh.triangle(id);
            triangles.push_back(Triangle(data.a, data.b, data.c, id));
        }

        printf("%d triangles\n", int(triangles.size()));
        assert(triangles.size());
    }

    Hit intersect(const Ray &ray) const
    {
        Hit hit;
        float tmax = ray.tmax;
        for (int id = 0; id < int(triangles.size()); id++)
            // ne renvoie vrai que si l'intersection existe dans l'intervalle [0 tmax]
            if (Hit h = triangles[id].intersect(ray, tmax))
            {
                hit = h;
                tmax = h.t;
            }

        return hit;
    }

    bool visible(const Ray &ray) const
    {
        for (int id = 0; id < int(triangles.size()); id++)
            if (triangles[id].intersect(ray, ray.tmax))
                return false;

        return true;
    }
};
*/
 
// struct BBoxNode{
//     Point pmin;
//     int left;
//     Point pmax;
//     int right;

//     BBoxNode BBoxNode { pmin, left,pmax, right };
// };



struct Node
{
    BBox bounds;
    int left;
    int right;
    int padleft;
    int padright;
    
    bool internal( ) const { return right > 0; }                        // renvoie vrai si le noeud est un noeud interne
    int internal_left( ) const { assert(internal()); return left; }     // renvoie le fils gauche du noeud interne 
    int internal_right( ) const { assert(internal()); return right; }   // renvoie le fils droit
    
    bool leaf( ) const { return right < 0; }                            // renvoie vrai si le noeud est une feuille
    int leaf_begin( ) const { assert(leaf()); return -left; }           // renvoie le premier objet de la feuille
    int leaf_end( ) const { assert(leaf()); return -right; }            // renvoie le dernier objet
};
 
// creation d'un noeud interne
Node make_node( const BBox& bounds, const int left, const int right )
{
    Node node { bounds, left, right };
    assert(node.internal());    // verifie que c'est bien un noeud...
    return node;
}
 
// creation d'une feuille
Node make_leaf( const BBox& bounds, const int begin, const int end )
{
    Node node { bounds, -begin, -end };
    assert(node.leaf());        // verifie que c'est bien une feuille...
    return node;
}



 struct triangle_less1
 {
     int axis;
     float cut;
     
     triangle_less1( const int _axis, const float _cut ) : axis(_axis), cut(_cut) {}
     
     bool operator() ( const Triangle& triangle ) const
     {
         // re-construit l'englobant du triangle
         BBox bounds= triangle.bounds();
         return bounds.centroid(axis) < cut;
     }
 };
  
  
 struct BVH
 {
     std::vector<Node> nodes;
     std::vector<Triangle> triangles;
     int root;
     
     int direct_tests;
     
     // construit un bvh pour l'ensemble de triangles
     int build( const BBox& _bounds, const std::vector<Triangle>& _triangles )
     {
         triangles= _triangles;  // copie les triangles pour les trier
         nodes.clear();          // efface les noeuds
         nodes.reserve(triangles.size());
         
         // construit l'arbre... 
         root= build(_bounds, 0, triangles.size());
         // et renvoie la racine
         return root;
     }
     
     void intersect( RayHit& ray ) const
     {
         Vector invd= Vector(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
         intersect( ray, invd);
     }
     
     void intersect_fast( RayHit& ray ) const
     {
         Vector invd= Vector(1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z);
         intersect_fast(root, ray, invd);
     }
     
 protected:
     // construction d'un noeud
     int build( const BBox& bounds, const int begin, const int end )
     {
         if(end - begin < 2)
         {
             // inserer une feuille et renvoyer son indice
             int index= nodes.size();
             nodes.push_back(make_leaf(bounds, begin, end));
             return index;
         }
         
         // axe le plus etire de l'englobant
         Vector d= Vector(bounds.pmin, bounds.pmax);
         int axis;
         if(d.x > d.y && d.x > d.z)  // x plus grand que y et z ?
             axis= 0;
         else if(d.y > d.z)          // y plus grand que z ? (et que x implicitement)
             axis= 1;
         else                        // x et y ne sont pas les plus grands...
             axis= 2;
  
         // coupe l'englobant au milieu
         float cut= bounds.centroid(axis);
         
         // repartit les triangles 
         Triangle *pm= std::partition(triangles.data() + begin, triangles.data() + end, triangle_less1(axis, cut));
         int m= std::distance(triangles.data(), pm);
         
         // la repartition des triangles peut echouer, et tous les triangles sont dans la meme partie... 
         // forcer quand meme un decoupage en 2 ensembles 
         if(m == begin || m == end)
             m= (begin + end) / 2;
         assert(m != begin);
         assert(m != end);
         
         // construire le fils gauche
         // les triangles se trouvent dans [begin .. m)
         BBox bounds_left= triangle_bounds(begin, m);
         int left= build(bounds_left, begin, m);
         
         // on recommence pour le fils droit
         // les triangles se trouvent dans [m .. end)
         BBox bounds_right= triangle_bounds(m, end);
         int right= build(bounds_right, m, end);
         
         int index= nodes.size();
         nodes.push_back(make_node(bounds, left, right));
         return index;
     }
     
     BBox triangle_bounds( const int begin, const int end )
     {
         BBox bbox= triangles[begin].bounds();
         for(int i= begin +1; i < end; i++)
             bbox.insert(triangles[i].bounds());
         
         return bbox;
     }
     
    //  void intersect( const int index, RayHit& ray, const Vector& invd ) const
    //  {
    //      const Node& node= nodes[index];
    //      if(node.bounds.intersect(ray, invd))
    //      {
    //          if(node.leaf())
    //          {
    //              for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
    //                  triangles[i].intersect(ray);
    //          }
    //          else // if(node.internal())
    //          {
    //              intersect(node.internal_left(), ray, invd);
    //              intersect(node.internal_right(), ray, invd);
    //          }
    //      }
    //  }

    void intersect( RayHit& ray, const Vector& invd ) const
     {

        int stack[64];
        int top= 0;
        
        // empiler la racine
        stack[top++]= root;
        
        //float tmax= ray.tmax;
        // tant qu'il y a un noeud dans la pile
        while(top > 0)
        {
            int index= stack[--top];
            
            const Node& node= nodes[index];
            if(node.bounds.intersect(ray, invd))
            {
                if(node.leaf())
                {
                    for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
                        triangles[i].intersect(ray);
                }
                else // if(node.internal())
                {
                    assert(top +1 < 64);       // le noeud est touche, empiler les fils
                    stack[top++]= node.internal_left();
                    stack[top++]= node.internal_right();
                }
            }
        }
     }
     void intersect_fast( const int index, RayHit& ray, const Vector& invd ) const
     {
         const Node& node= nodes[index];
         if(node.leaf())
         {
             for(int i= node.leaf_begin(); i < node.leaf_end(); i++)
                 triangles[i].intersect(ray);
         }
         else // if(node.internal())
         {
             const Node& left_node= nodes[node.left];
             const Node& right_node= nodes[node.right];
             
             BBoxHit left= left_node.bounds.intersect(ray, invd);
             BBoxHit right= right_node.bounds.intersect(ray, invd);
             if(left && right)                                                   // les 2 fils sont touches par le rayon...
             {
                 if(left.centroid() < right.centroid())                          // parcours de gauche a droite
                 {
                     intersect_fast(node.internal_left(), ray, invd);
                     intersect_fast(node.internal_right(), ray, invd);
                 }
                 else                                                            // parcours de droite a gauche                                        
                 {
                     intersect_fast(node.internal_right(), ray, invd);
                     intersect_fast(node.internal_left(), ray, invd);
                 }
             }
             else if(left)                                                       // uniquement le fils gauche
                 intersect_fast(node.internal_left(), ray, invd);
             else if(right)
                 intersect_fast(node.internal_right(), ray, invd);               // uniquement le fils droit
         }
     } 
 };
  

struct NodeGPU
{
    vec3 pmin;
    int left;
    vec3 pmax;
    int right;

    NodeGPU( const Node& node) : pmin(node.bounds.pmin),left(node.left), pmax(node.bounds.pmax),right(node.right) {}
    
};



struct RT : public AppTime
{
    // constructeur : donner les dimensions de l'image, et eventuellement la version d'openGL.
    RT( const char *filename ) : AppTime(1024, 640) 
    {
        m_mesh= read_mesh(filename);
    }
    
    int init( )
    {
        if(m_mesh.triangle_count() == 0)
            return -1;
        
        Point pmin, pmax;
        m_mesh.bounds(pmin, pmax);
        m_camera.lookat(pmin, pmax);
        
        // utilise un compute shader qui :
        // 1. cree un rayon pour chaque pixel de l'image
        // 2. calcule les intersections du rayon avec tous les triangles. (et garde la plus proche...)
        // 3. ecrit le resultat dans une image.
        
        // il faut donc creer le buffer et la texture, les entrees / sorties du shader
        
        // recupere les triangles du mesh
        // structure declaree par le shader, en respectant l'alignement std430


        BBox bounds;
        m_mesh.bounds(bounds.pmin, bounds.pmax);
        
        
        std::vector<Triangle> data;
        data.reserve(m_mesh.triangle_count());
        for(int i= 0; i < m_mesh.triangle_count(); i++)
        {
            TriangleData t= m_mesh.triangle(i);
            data.push_back( { Point(t.a),1, Point(t.b) ,1, Point(t.c) , i } );
            //Triangle( const Point &_a, const int _padp, const Point &_b, const int _pade1 , const Point &_c, const int _id) 
            //data.emplace_back(t, i);
        }
 
        BVH bvh;
        
        {
            auto start= std::chrono::high_resolution_clock::now();
            // construction 
            bvh.build(bounds, data);
            
            auto stop= std::chrono::high_resolution_clock::now();
            int cpu= std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
            printf("build %dms\n", cpu);
        }

        root_uni = bvh.root;
        std::cout<<"root_uni "<<root_uni<<std::endl;


        std::vector<Node> dataNode;
        //data.reserve(m_mesh.triangle_count());
        for(int i= 0; i < bvh.nodes.size() ; i++)
        {
            //Node BNode= Node(bvh.nodes[i].bounds, bvh.nodes[i].left, bvh.nodes[i].right, bvh.nodes[i].padleft, bvh.nodes[i].padright);
            Node BNode= Node(bvh.nodes[i]);
            //data.push_back( { Point(t.a),1, Point(t.b) - Point(t.a),1, Point(t.c) - Point(t.a), i } );
            dataNode.emplace_back(BNode);
        }

        std::vector<NodeGPU> dataNodeGPU;
        //data.reserve(m_mesh.triangle_count());
        for(int i= 0; i < bvh.nodes.size(); i++)
        {
            NodeGPU BNode= NodeGPU(bvh.nodes[i]);
            //data.push_back( { Point(t.a),1, Point(t.b) - Point(t.a),1, Point(t.c) - Point(t.a), i } );
            dataNodeGPU.emplace_back(BNode);
        }

        std::vector<Triangle> trianglesGPU;
        {
            int n= bvh.triangles.size();
            std::cout<<"n : "<<n<<std::endl;
            for(int i= 0; i < n; i++){
                trianglesGPU.emplace_back(bvh.triangles[i]);
                //std::cout<<"triangles id GPU: "<<trianglesGPU[i].id<<std::endl;

            }
        }

        std::vector<Triangle> trianglesGPU2;
        data.reserve(m_mesh.triangle_count());
        for(int i= 0; i < m_mesh.triangle_count(); i++)
        {
            TriangleData t= m_mesh.triangle(i);
            trianglesGPU2.push_back( { Point(t.a),1, Point(t.b) ,1, Point(t.c) , bvh.triangles[i].id } );
            //std::cout<<"triangles id GPU2: "<<bvh.triangles[i].id<<std::endl;

            //Triangle( const Point &_a, const int _padp, const Point &_b, const int _pade1 , const Point &_c, const int _id) 
            //data.emplace_back(t, i);
        }

        // std::vector<Node> nodes;
        // data.reserve(m_mesh.triangle_count());
        // for(int i= 0; i < m_mesh.triangle_count(); i++)
        // {
        //     TriangleData t= m_mesh.triangle(i);
        //     data.push_back( { Point(t.a), Point(t.b) - Point(t.a), Point(t.c) - Point(t.a) } );
        // }
        
        // cree et initialise le storage buffer 0
        glGenBuffers(1, &m_buffer_tri);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_buffer_tri);
        glBufferData(GL_SHADER_STORAGE_BUFFER, trianglesGPU.size() * sizeof(Triangle), trianglesGPU.data(), GL_STATIC_READ);

        // // cree et initialise le storage buffer 0
        // glGenBuffers(1, &m_buffer_tri);
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_buffer_tri);
        // glBufferData(GL_SHADER_STORAGE_BUFFER, data.size() * sizeof(Triangle), data.data(), GL_STATIC_READ);
        
        // cree et initialise le storage buffer 1
        glGenBuffers(1, &m_buffer_node);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_buffer_node);
        glBufferData(GL_SHADER_STORAGE_BUFFER, bvh.nodes.size() * sizeof(Node), dataNode.data(), GL_STATIC_READ);

        // // cree et initialise le storage buffer 1
        // glGenBuffers(1, &m_buffer_node);
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_buffer_node);
        // glBufferData(GL_SHADER_STORAGE_BUFFER, dataNodeGPU.size() * sizeof(NodeGPU), dataNodeGPU.data(), GL_STATIC_READ);
    
        // // cree et initialise le storage buffer 2
        // glGenBuffers(1, &m_buffer_tri2);
        // glBindBuffer(GL_SHADER_STORAGE_BUFFER, m_buffer_tri2);
        // glBufferData(GL_SHADER_STORAGE_BUFFER, data.size() * sizeof(Triangle), data.data(), GL_STATIC_READ);
        

        // texture / image resultat
        // cree la texture, 4 canaux, entiers 8bits normalises, standard
        glGenTextures(1, &m_texture);
        glBindTexture(GL_TEXTURE_2D, m_texture);
        glTexImage2D(GL_TEXTURE_2D, 0,
            GL_RGBA8, window_width(), window_height(), 0,
            GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        
        // pas la peine de construire les mipmaps, le shader ne va ecrire que le mipmap 0
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
        // oui, c'est une texture tout a fait normale.




        int wh = window_height()*window_width();
        std::random_device seed_generator;
        std::vector<unsigned int> seeds(wh);
        for (int i = 0; i< wh; i++ ){
            seeds[i] = seed_generator();
        }



        // texture / image resultat
        // cree la texture, 4 canaux, entiers 8bits normalises, standard
        glGenTextures(1, &m_seed_image);
        glBindTexture(GL_TEXTURE_2D, m_seed_image);
        glTexImage2D(GL_TEXTURE_2D, 0,
            GL_R32UI, window_width(), window_height(), 0,
            GL_RED_INTEGER, GL_UNSIGNED_INT, seeds.data());
        
        // pas la peine de construire les mipmaps, le shader ne va ecrire que le mipmap 0
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);
        // oui, c'est une texture tout a fait normale.

        // pour afficher l'image resultat, 2 solutions : 
        //      1. utiliser un shader qui copie un pixel de la texture vers un pixel du framebuffer par defaut / la fenetre, 
        //      2. ou copier directement la texture sur le framebuffer par defaut / la fenetre, en utilisant glBlitFramebuffer
        
        // pour changer, on va utiliser glBlitFramebuffer, 
        // mais il faut configurer un framebuffer...
        glGenFramebuffers(1, &m_blit_framebuffer);
        glBindFramebuffer(GL_READ_FRAMEBUFFER, m_blit_framebuffer);
        glFramebufferTexture(GL_READ_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, m_texture, 0);
        // selectionner la texture du framebuffer a copier...
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        
        // compute shader
        m_program= read_program("tutos/M2/raytrace_compute.glsl");
        program_print_errors(m_program);

        frame=0;
        
        // nettoyage
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindTexture(GL_TEXTURE_2D, 0);
        return 0;
    }
    
    int quit( )
    {
        release_program(m_program);
        glDeleteTextures(1, &m_texture);
        glDeleteTextures(1, &m_seed_image);
        glDeleteBuffers(1, &m_buffer_tri);
        glDeleteBuffers(1, &m_buffer_node);
        glDeleteFramebuffers(1, &m_blit_framebuffer);
        return 0;
    }
    
    int render( )
    {
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glViewport(0, 0, window_width(), window_height());
        glClear(GL_COLOR_BUFFER_BIT);
        
        if(key_state('f'))
        {
            clear_key_state('f');
            // recentrer la camera
            Point pmin, pmax;
            m_mesh.bounds(pmin, pmax);
            m_camera.lookat(pmin, pmax);        
        }
        
        // deplace la camera
        int mx, my;
        unsigned int mb= SDL_GetRelativeMouseState(&mx, &my);
        if(mb & SDL_BUTTON(1)){           // le bouton gauche est enfonce
            m_camera.rotation(mx, my);
            frame=0;
        } 
        else if(mb & SDL_BUTTON(3)){         // le bouton droit est enfonce
            m_camera.move(mx); 
            frame=0;
        } 
        else if(mb & SDL_BUTTON(2)){     // le bouton du milieu est enfonce
            m_camera.translation((float) mx / (float) window_width(), (float) my / (float) window_height());
            frame=0;
        }    
        
        
        // recupere les transformations standards.
        Transform m= Identity();
        Transform v= m_camera.view();
        Transform p= m_camera.projection(window_width(), window_height(), 45);
        Transform im= Viewport(window_width(), window_height());
        // ou Transform im= m_camera.viewport();
        // compose toutes les transformations, jusqu'au repere image
        Transform T= im * p * v * m;
        
        // config pipeline
        glUseProgram(m_program);
        
        // storage buffer 0
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, m_buffer_tri);

        // storage buffer 1
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, m_buffer_node);

        // // storage buffer 2
        // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, m_buffer_tri2);

        
        // image texture 0, ecriture seule, mipmap 0 + format rgba8 classique
        glBindImageTexture(0, m_texture, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA8);
        // configurer le shader
        program_uniform(m_program, "image", 0);

        glBindImageTexture(1, m_seed_image, 0, GL_FALSE, 0, GL_READ_ONLY, GL_R32UI);
        // configurer le shader
        program_uniform(m_program, "seeds", 1);
        
        // uniforms
        program_uniform(m_program, "invMatrix", T.inverse());

        program_uniform(m_program, "root", root_uni);

        program_uniform(m_program, "frame", frame);
        //glUniform1i( location(m_program, "frame"), frame)
        frame++;
        
        // nombre de groupes de shaders pour executer un compute shader par pixel de l'image resultat. on utilise un domaine 2d...
        // le shader declare un groupe de threads de 8x8.
        int nx= window_width() / 8;
        int ny= window_height() / 8;
        // on suppose que les dimensions de l'image sont multiples de 8...
        // sinon calculer correctement le nombre de groupes pour x et y. 
        
        // go !!
        glDispatchCompute(nx, ny, 1);
        
        // attendre le resultat
        glMemoryBarrier(GL_ALL_BARRIER_BITS);
        
        // afficher le resultat
        // copier la texture resultat vers le framebuffer par defaut / de la fenetre
        glBindFramebuffer(GL_READ_FRAMEBUFFER, m_blit_framebuffer);
        
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        glBlitFramebuffer(
            0,0, window_width(),window_height(), 
            0,0, window_width(),window_height(),
             GL_COLOR_BUFFER_BIT, GL_NEAREST);
        
        return 1;
    }
    
protected:
    Mesh m_mesh;
    Orbiter m_camera;

    GLuint m_blit_framebuffer;
    GLuint m_program;
    GLuint m_texture;
    GLuint m_seed_image;
    GLuint m_buffer_tri;
    GLuint m_buffer_node;

    int frame;
    int root_uni;
};

    
int main( int argc, char **argv )
{
    const char *filename= "data/cornell.obj";
    //const char *filename= "data/247_House_15_obj.obj";
    //const char *filename= "data/emission.obj";
    if(argc > 1)
        filename= argv[1];
    
    RT app(filename);
    app.run();
    
    return 0;
}
