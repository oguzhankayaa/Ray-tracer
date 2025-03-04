#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <thread>
//#include <chrono>

class ray
{
public:
    parser::Vec3f o, d;
    int depth=0;
};

struct hitRecord
{
    int ObjI=-1;
    parser::Vec3f hitPoint,d,normal;
    float t=-1;
    parser::Material mat;
    bool triangle= false;
    
};




parser::Vec3f crossProduct(const parser::Vec3f& a, const parser::Vec3f& b)
{
    parser::Vec3f result;
    result.x = a.y*b.z-a.z*b.y;
    result.y = a.z*b.x-a.x*b.z;
    result.z = a.x*b.y-a.y*b.x;
    return result;
}


float dotProduct(const parser::Vec3f& a, const parser::Vec3f& b)
{
    return a.x*b.x+ a.y*b.y+ a.z*b.z;
}


parser::Vec3f multS(const parser::Vec3f& a, float b)
{
    parser::Vec3f result;
    result.x = a.x*b;
    result.y = a.y*b;
    result.z = a.z*b;
    return result;
}


parser::Vec3f addS(const parser::Vec3f& a, float b)
{
    parser::Vec3f result;
    result.x = a.x+b;
    result.y = a.y+b;
    result.z = a.z+b;
    return result;
}



parser::Vec3f add(const parser::Vec3f& a, const parser::Vec3f& b)
{
    parser::Vec3f result;
    result.x = a.x+b.x;
    result.y = a.y+b.y;
    result.z = a.z+b.z;
    return result;
}



parser::Vec3f weirdProduct(const parser::Vec3f& a, const parser::Vec3f& b)
{
    parser::Vec3f result;
    result.x = a.x*b.x;
    result.y = a.y*b.y;
    result.z = a.z*b.z;
    return result;
}


parser::Vec3f normalize(parser::Vec3f a)
{
	return multS(a,1.0/sqrt(dotProduct(a,a)));
}


parser::Vec3f truncate(parser::Vec3f a, float min, float max)
{
	if (a.x>max)
    {
        a.x=max;
    }
    else if (a.x<min)
    {
        a.x=min;
    }
    
    if (a.y>max)
    {
        a.y=max;
    }
    else if (a.y<min)
    {
        a.y=min;
    }
    if (a.z>max)
    {
        a.z=max;
    }
    else if (a.z<min)
    {
        a.z=min;
    }
    return a;
    
}

float intersectSphere(ray r, parser::Sphere s ,parser::Scene &scene)
{
    float A,B,C,delta;
    parser::Vec3f c;
    c=scene.vertex_data[s.center_vertex_id-1];

    float t1,t2,t;
    
    C=(dotProduct(add(r.o,multS(c,-1)),add(r.o,multS(c,-1))))-(s.radius*s.radius);
    B=2*dotProduct(r.d,add(r.o,multS(c,-1)));
    A=dotProduct(r.d,r.d);
    delta=(B*B)-(4*A*C);
    //std::cout<<"A"<<A<<std::endl;
    //std::cout<<"B"<<B<<std::endl;
    //std::cout<<"C"<<C<<std::endl;

    //std::cout<<"delta:"<<delta<<std::endl;
    if (delta<0) return -1;
    else if(delta==0)
    {
        t=-B/(2*A);
        //std::cout<<"0"<<std::endl;
    }
    
    else
    {
       //std::cout<<"delta büyük"<<std::endl;
        t1=(-B+ sqrt(delta))/(2*A);
        t2=(-B- sqrt(delta))/(2*A);
        t=t1;
        if (t1>t2)
        {
            t=t2;
        }
        
    }
    //std::cout<<"t: "<<t<<std::endl;
    return t;


}

float intersectTriangle(ray r, parser::Face face, parser::Scene &scene)
{
    parser::Vec3f a,b,c,a_b,a_c,a_o;
    float beta, gamma, t,deta;
    a=scene.vertex_data[face.v0_id-1];
    b=scene.vertex_data[face.v1_id-1];
    c=scene.vertex_data[face.v2_id-1];

    a_b= add(a,multS(b,-1));  
    a_c= add(a,multS(c,-1));  
    a_o= add(a,multS(r.o,-1));

    deta=dotProduct(a_b,crossProduct(a_c,r.d));
        
    
    t=dotProduct(a_b,crossProduct(a_c,a_o)) /deta;
    gamma=dotProduct(a_b,crossProduct(a_o,r.d)) /deta;
    beta=dotProduct(a_o,crossProduct(a_c,r.d)) /deta;

    if (beta >= 0 && gamma >= 0 && beta + gamma <= 1 && t > 0)
    {
        return t;
    }
    return -1;
}
/////****************End of helpers */

ray generateRay(int j , int i, parser::Camera cam)
{

    float left = cam.near_plane.x;
    float right = cam.near_plane.y;
    float bottom = cam.near_plane.z;
    float top = cam.near_plane.w;
    ray result;
    parser::Vec3f w, u, m, q,s;
    float su, sv;
    su = (i+0.5)*(right-left)/cam.image_width;
    sv = (j+0.5)*(top-bottom)/cam.image_height; 
    w=multS(cam.gaze, -1);
    u=crossProduct(cam.up,w);
    m=add(cam.position,multS(cam.gaze,cam.near_distance));//negative?
    q=add(m,add(multS(u,left),multS(cam.up,top)));
    s=add(q,add(multS(u,su),multS(cam.up,-sv)));
    
    result.o= cam.position;
    result.d= add(s,multS(cam.position,-1));
    result.d=normalize(result.d);///should I normalize that?
    return result;
}

bool closestHit(ray r, hitRecord& hr, parser::Scene &scene)
{ 
    parser::Sphere sph;
    parser::Mesh msh;
    parser::Triangle tr;
    float t, mshminT=999999.000,trminT=999999.000,sphminT=999999.000;
    int mshminI=-1,trminI=-1, sphminI=-1;
    for (int i = 0; i < scene.spheres.size() ; i++)
    {
        t=intersectSphere(r, scene.spheres[i],scene);
        //std::cout<<"t: "<<r.depth<<std::endl;
        if (t<sphminT && t>=0)//what should be t greater than?
        {   
            sphminT=t;
            sphminI=i;
           //std::cout<<"t: "<<t<<std::endl;
           //std::cout<<"i: "<<i<<std::endl;
           //std::cout<<"matid: "<<colorMatId<<std::endl;
        }
        
    }
    for (int i = 0; i < scene.triangles.size() ; i++)
    {
        t=intersectTriangle(r, scene.triangles[i].indices,scene);
        //std::cout<<"t: "<<r.depth<<std::endl;
        if (t<trminT && t>=0)//what should be t greater than?
        {   
            trminT=t;
            trminI=i;
           //std::cout<<"t: "<<t<<std::endl;
           //std::cout<<"i: "<<i<<std::endl;
           //std::cout<<"matid: "<<colorMatId<<std::endl;
        }
        
    }
    
    
    for (int i = 0; i < scene.meshes.size(); i++)
    {  
        for (int k = 0; k < scene.meshes[i].faces.size(); k++)
        {
            t=intersectTriangle(r, scene.meshes[i].faces[k],scene);
            if (t<trminT && t>=0)//what should be t greater than?
            {   
                trminT=t;
                trminI=k;
                mshminI=i;
               //std::cout<<"t: "<<t<<std::endl;
               //std::cout<<"i: "<<i<<std::endl;
               //std::cout<<"matid: "<<colorMatId<<std::endl;
            }
        }
        
    }
    
    
    if (trminT<sphminT)
    {
        parser::Face face;
        if (mshminI==-1)
        {
            hr.mat=scene.materials[scene.triangles[trminI].material_id-1];
            face=scene.triangles[trminI].indices;
        }
        else
        {
            hr.mat=scene.materials[scene.meshes[mshminI].material_id-1];
            face=scene.meshes[mshminI].faces[trminI];
        }
        
        hr.hitPoint = add(r.o,multS(r.d,trminT));
        hr.ObjI = trminI;
        hr.triangle=true;
        hr.t = trminT;
        hr.d = add(r.o,multS(hr.hitPoint,-1));//should I normalize
        /*
        parser::Vec3f a,b,c,s1,s2;
        a=scene.vertex_data[face.v0_id-1];
        b=scene.vertex_data[face.v1_id-1];
        c=scene.vertex_data[face.v2_id-1];
        s1=add(b,multS(a,-1)); 
        s2=add(c,multS(a,-1)); 


        hr.normal=crossProduct(s1,s2);//center ID 0 dan mı 1 den mi başlar
        hr.normal=normalize(hr.normal);*/
        hr.normal=face.normal;
        return true;
    }
    

    if (sphminI!=-1)
    {
        hr.hitPoint = add(r.o,multS(r.d,sphminT));
        hr.ObjI = sphminI;
        hr.mat=scene.materials[scene.spheres[sphminI].material_id-1];
        hr.t = sphminT;
        hr.d = add(r.o,multS(hr.hitPoint,-1));//should I normalize

        int centerID=scene.spheres[sphminI].center_vertex_id;
        hr.normal=add(hr.hitPoint,multS(scene.vertex_data[centerID-1],-1));//center ID 0 dan mı 1 den mi başlar
        hr.normal=normalize(hr.normal);

        /*if (r.depth>0)
        {
            printf("%d----------------------------------------------\n",r.depth);
        }*/
        
        return true;
    }
    //printf("%d\n",r.depth);
    return false;
    
    
}

bool inShadow(hitRecord hr, parser::PointLight light, parser::Scene &scene)
{
    float t, distance;
    bool shadow=false;
    parser::Vec3f wi,L;
    wi=add(light.position ,multS(hr.hitPoint,-1));//ray in. from hit to light
    
    L=normalize(wi);
    ray lightRay;
    lightRay.o=hr.hitPoint;
    lightRay.o=add(lightRay.o,multS(hr.normal, scene.shadow_ray_epsilon));///TODO : what should be the direction? maybe L
    lightRay.d=L;
    distance=sqrt(dotProduct(wi, wi));
    hitRecord hr2;
    closestHit(lightRay,hr2,scene);
    if (hr2.t<distance && hr2.t>0)//what should be t greater than?
    {   
        shadow=true;
    }
    return shadow;

}

parser::Vec3f diffuse(hitRecord hr, parser::PointLight light)
{
    parser::Vec3f wi,L,diffuse;
    diffuse=hr.mat.diffuse;
    wi=add(light.position ,multS(hr.hitPoint,-1));//ray out
    L=normalize(wi);
    float dot=dotProduct(hr.normal,L);
    if (dot>0)// TODO do we need that? 
    {
        diffuse=multS(diffuse,dot);
        //std::cout<<"diff1: "<<diffuse.x<<" "<<diffuse.y<<" "<<diffuse.z<<" "<<std::endl;
        diffuse=weirdProduct(diffuse,light.intensity);
        //std::cout << "LightDistance: " << dotProduct(wi,wi) << std::endl;
        diffuse=multS(diffuse,1/(dotProduct(wi,wi)));
    }
    //else std::cout << "bu " << std::endl;
    else{
        diffuse={0,0,0};
    }
    return diffuse;
}

parser::Vec3f spec(hitRecord hr, parser::PointLight light)
{   
    parser::Vec3f wi,L,specular,H,color;
    float cosalph;
    specular=hr.mat.specular;
    wi=add(light.position ,multS(hr.hitPoint,-1));//ray out
    L=normalize(wi);
    float dot=dotProduct(hr.normal,L);
    color.x=color.y=color.z=0;
    if (dot>0)
    {
        H=add(L,normalize(hr.d));//wi+wo  TODO  should it be -1 or direct addition
        H=normalize(H);
        cosalph=dotProduct(H,hr.normal);
        if (cosalph<0)
        {
            cosalph=0;
        }
        color=weirdProduct(multS(light.intensity,pow(cosalph,hr.mat.phong_exponent)/(dotProduct(wi,wi))),specular);
        
    }
    return color;

}
parser::Vec3f computeColor(ray r, parser::Scene &scene);

parser::Vec3f applyShading(ray r, hitRecord hr,parser::Scene &scene)
{
    parser::Vec3f color= weirdProduct(hr.mat.ambient,scene.ambient_light);
    if (hr.mat.is_mirror)
    {   
        ray reflectionRay;
        parser::Vec3f Wr;
        Wr=add(multS(normalize(hr.d),-1),multS(normalize(hr.normal),2*dotProduct(normalize(hr.normal),normalize(hr.d))));
        Wr=normalize(Wr);
        reflectionRay.depth=r.depth+1;
        //std::cout<<r.depth<<std::endl;
        reflectionRay.o=hr.hitPoint;
        reflectionRay.o = add(hr.hitPoint, multS(hr.normal, scene.shadow_ray_epsilon));
        reflectionRay.d=Wr;
        color=add(color,weirdProduct(computeColor(reflectionRay,scene),hr.mat.mirror));
        /* inşallah yazarız */
    }

    for ( int i = 0; i < scene.point_lights.size(); i++)
    {
        if (! inShadow(hr,scene.point_lights[i],scene))
        {
            color=add(color,diffuse(hr,scene.point_lights[i]));
            color=add(color,spec(hr,scene.point_lights[i]));
            
        }
        
    }

    return color;
}


parser::Vec3f computeColor(ray r,parser::Scene &scene)
{
    parser::Vec3f bc;
    bc.x=scene.background_color.x;
    bc.y=scene.background_color.y;
    bc.z=scene.background_color.z;
    
    parser::Vec3f zero;
    zero.x=zero.y=zero.z=0;
    if (r.depth>scene.max_recursion_depth)
    {
        
        return zero;
    }
    hitRecord hr;
    if (closestHit(r,hr,scene))
    {
       // std::cout<<r.depth<<std::endl;
        return applyShading(r,hr,scene);
    }
    else if (r.depth==0)
    {
        return bc;
    }
    else
    {
        return zero;
    }
}
void division(int start, int end, unsigned char* image, parser::Scene scene, parser::Camera cam)
{
    int Pixis=3*start*cam.image_width;
     for (int y = start; y < end; ++y)
        {
            for (int x = 0; x < cam.image_width; ++x)
            {
                ray myRay= generateRay(y,x,cam);
                //std::cout<<"ox:"<<myRay.o.x<<"oy:"<<myRay.o.y<<std::endl;  
                //std::cout<<"dx:"<<myRay.d.x<<"dy:"<<myRay.d.y<<std::endl;
                parser::Vec3f rayC;
                rayC=computeColor(myRay,scene);
                rayC = truncate(rayC, 0.0f, 255.0f);
                image[Pixis++] = rayC.x;
                image[Pixis++] = rayC.y;
                image[Pixis++] = rayC.z;
               /* if (rayC.y!=0)
                {
                    std::cout<<y<<" "<<x<<" green: "<<rayC.y<<std::endl;
                }*/
                
                

            //std::cout<<rayC.x<<std::endl;
            }
        }
}

void computeNormal(parser::Face &face, parser::Scene &scene)
{
    parser::Vec3f a,b,c,s1,s2;
    a=scene.vertex_data[face.v0_id-1];
    b=scene.vertex_data[face.v1_id-1];
    c=scene.vertex_data[face.v2_id-1];
    s1=add(b,multS(a,-1)); 
    s2=add(c,multS(a,-1)); 


    face.normal=normalize((crossProduct(s1,s2)));//center ID 0 dan mı 1 den mi başlar

}

void trianglePre(parser::Scene &scene)
{   
    for (int i = 0; i < scene.meshes.size(); i++)
    {
        for (int k= 0;  k< scene.meshes[i].faces.size(); k++)
        {
            computeNormal(scene.meshes[i].faces[k],scene);  
        }
        
    }
    for (int i = 0; i < scene.triangles.size(); i++)
    {
        computeNormal(scene.triangles[i].indices,scene);
    }
    return;
    

}
int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    trianglePre(scene);
    parser::Vec3i backgroundColor = scene.background_color;

    //auto start = std::chrono::high_resolution_clock::now();

    for (int camNum = 0; camNum < scene.cameras.size(); camNum++)
    {
        parser::Camera cam= scene.cameras[camNum];
        int width = cam.image_width, height =cam.image_height;

        unsigned char* image = new unsigned char [width * height * 3];

        //setBackground(image,backgroundColor, width, height);
        //main ray tracing

int end1i = cam.image_height /4;
int end2i = end1i+ cam.image_height / 4;
int end3i = end2i+ cam.image_height / 4;
int end4i = cam.image_height;
std::thread t1(division, 0, end1i, image, scene, cam);
std::thread t2(division, end1i, end2i, image, scene, cam);
std::thread t3(division, end2i, end3i, image, scene, cam);
std::thread t4(division, end3i, end4i, image, scene, cam);
t1.join();
t2.join();
t3.join();
t4.join();


        
            std::cout<<cam.image_name<<std::endl;
            write_ppm(&cam.image_name[0], image, width, height);

    }

    /*auto end = std::chrono::high_resolution_clock::now();

    // Calculate the elapsed time in minutes
    std::chrono::duration<double> elapsed = end - start;
    double elapsedMinutes = elapsed.count() / 60.0;

    // Output the elapsed time
    std::cout << "Time taken: " << elapsedMinutes << " minutes" << std::endl;*/

    return 0;
}