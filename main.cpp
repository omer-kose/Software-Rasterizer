#include <vector>
#include <iostream>
#include <algorithm>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "GeneralGraphicsUtilites.h"

Model *model = nullptr;
float *shadowbuffer = nullptr;

const int width = 800;
const int height = 800;

/*Actually this one is in the reverse direction of the direction, we intended*/
Vec3f light_dir(1, 1, 0);


Vec3f       eye(1, 1, 4);
Vec3f    center(0, 0, 0);
Vec3f        up(0, 1, 0);

/*Model Properties*/
Vec3f scalingFactors(1.,1.,1.0);
Vec3f position(0.0, 0.0, 0.0);
double pitch = 0.0;
double yaw = -60.0;
double roll =0.0;

/*External Geometry.h notes*/
/*Embed<DIM>*/
/*embed function takes a vector of dimension n and returns another vector of dimension templated*/
/*embed<4> returns a 4x1 column vector for example, It fills rest of the columns with a second parameter (defaulted to 1)*/
/*It is used to transform vectors to the homogeneous coordinates*/


/*Proj<dim>*/
/*Similiar to embed but instead it lowers the dimension of the vector*/
/*It is used to transform from homogeneous to cartesian*/


struct DepthShader : public IShader
{
	mat<3, 3, float> varying_tri;

	DepthShader () : varying_tri() {}

	virtual Vec4f vertex(int iface, int nthvert)
	{
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_Vertex = Viewport * Projection*View * World *gl_Vertex;          // transform it to screen coordinates
		varying_tri.set_col(nthvert, proj<3>(gl_Vertex / gl_Vertex[3]));
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) 
	{
		Vec3f p = varying_tri * bar;
		color = TGAColor(255, 255, 255)*(p.z / depth);
		return false;
	}

};


struct SpecularWithNormalMappingShader : public IShader {
	mat<2, 3, float> varying_uv;  // same as above
	mat<4, 4, float> uniform_M;   //  World (might also need to include View not sure for now)
	mat<4, 4, float> uniform_MIT; // (World).invert_transpose() (might also need to include View not sure for now)

	virtual Vec4f vertex(int iface, int nthvert) {
		varying_uv.set_col(nthvert, model->uv(iface, nthvert));
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		return Viewport * Projection*View * World * gl_Vertex; // transform it to screen coordinates
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) 
	{
		Vec2f uv = varying_uv * bar;
		Vec3f normal = proj<3>(uniform_MIT*embed<4>(model->normal(uv))).normalize();
		Vec3f reflected = ((normal * (2 * (normal * light_dir))) - light_dir).normalize();
		float spec = std::max(0.f, std::pow(eye.normalize() * reflected, 1));
		float diff = std::max(0.f, normal*light_dir);
		TGAColor c = model->diffuse(uv);
		color = c;
		for (int i = 0; i < 3; i++) color[i] = std::min<float>(0.35 + c[i] * (diff + 0.5*spec), 255);
		return false;
	}
};


struct NormalMappingShader : public IShader 
{
	mat<2, 3, float> varying_uv;  // same as above
	/*Uniform mimics the OpenGL's uniform keyboard*/
	/*Uniform variables allows us to pass constants from outside to the shaders*/
	mat<4, 4, float> uniform_M;   //  Projection* View * World
	mat<4, 4, float> uniform_MIT; // (M).invert_transpose()

	virtual Vec4f vertex(int iface, int nthvert)
	{
		varying_uv.set_col(nthvert, model->uv(iface, nthvert));
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		return Viewport * Projection* View * World *gl_Vertex; // transform it to screen coordinates
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) 
	{
		Vec2f uv = varying_uv * bar;                 // interpolate uv for the current pixel
		/*Fill normals from normal map*/
		/*Normals are stored wrt local space we should take inverse transpose of the transformation(composition of matrices actually) matrix*/
		/*That we applied to the corresponding point*/
		Vec3f n = proj<3>(uniform_MIT*embed<4>(model->normal(uv))).normalize();
		float intensity = std::max(0.f, n*light_dir);
		/*Colors are stored in diffuse textures(maps)*/
		color = model->diffuse(uv)*intensity;      
		return false;                              
	}
};

/*Here texture is actually a diffuse texture*/
struct TextureShader : public IShader
{
	Vec3f varying_intensity;/*Mimic to OpenGL's reserved keyword varying. Written by vertex shader, read by fragment shader*/
	mat<2, 3, float> varying_uv;/*Same logic applies here, use barycentric coordinates to interpolate texture colors*/

	/*Iface is the triangle index, nthvert is the nth vertex of that triangle*/
	virtual Vec4f vertex(int iface, int nthvert)
	{
		varying_uv.set_col(nthvert, model->uv(iface, nthvert));
		Vec3f normal = model->normal(iface, nthvert);
		Vec3f reflected = (normal * (2 * (normal * light_dir))) - light_dir;
		varying_intensity[nthvert] = std::max(0.0f, std::pow(eye.normalize() * reflected, 1));
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));
		return Viewport * Projection *View * World * gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) 
	{
		float intensity = (varying_intensity * bar) + 0.35;   // interpolate intensity for the current pixel
		Vec2f uv = varying_uv * bar;                 // interpolate uv for the current pixel
		color = model->diffuse(uv)*intensity;      
		return false;                              
	}
};



/*Useless*/
struct DielectricShader : public IShader
{
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	virtual Vec4f vertex(int iface, int nthvert) {
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_Vertex = Viewport * Projection * View * World *gl_Vertex;     // transform it to screen coordinates
		/*This data is sent to fragment shader*/
		/*Compute the reflected vector*/
		Vec3f normal = model->normal(iface, nthvert);
		double etai = 1.0; /*Air*/
		double etat = 2.0; /*Some material*/
		/*Eta i over eta t*/
		double eiot = etai / etat;
		double cosI = light_dir * normal;
		double sinO = sqrt(1 - (eiot * eiot) * (1 - (cosI * cosI)));
		Vec3f refracted = (((normal * cosI) - light_dir) - normal * sinO) * eiot;
		varying_intensity[nthvert] = std::max(0.f, normal * (-1) * refracted);
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		float intensity = varying_intensity * bar;   // interpolate intensity for the current pixel
		color = TGAColor(255, 255, 255)*intensity; // well duh
		return false;                              // no, we do not discard this pixel
	}
};


struct PointLightShader : public IShader
{
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	virtual Vec4f vertex(int iface, int nthvert) {
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_Vertex = Viewport * Projection * View * World *gl_Vertex;     // transform it to screen coordinates
		/*This data is sent to fragment shader*/
		/*Compute the reflected vector*/
		Vec3f normal = model->normal(iface, nthvert);
		Vec3f reflected = (normal * (2 * (normal * light_dir))) - light_dir;
		/*Using light_dir as light position*/
		double attenuation = 1.0 / ((model->vert(iface, nthvert) - light_dir).norm() * (model->vert(iface, nthvert) - light_dir).norm());
		varying_intensity[nthvert] = std::max(0.0, std::pow(eye.normalize() * reflected * attenuation, 1));
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		float intensity = (varying_intensity * bar);   // interpolate intensity for the current pixel
		color = TGAColor(255, 255, 255)*intensity; // well duh
		return false;                              // no, we do not discard this pixel
	}
};



struct DirectionalLightShader : public IShader
{
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	virtual Vec4f vertex(int iface, int nthvert) {
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_Vertex = Viewport * Projection * View * World *gl_Vertex;     // transform it to screen coordinates
		/*This data is sent to fragment shader*/
		/*Compute the reflected vector*/
		Vec3f normal = model->normal(iface, nthvert);
		Vec3f reflected = ((normal * (2 * (normal * light_dir)))- light_dir).normalize();
		varying_intensity[nthvert] = std::max(0.f, std::pow(eye.normalize() * reflected,1));
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		float intensity = varying_intensity * bar;   // interpolate intensity for the current pixel
		color = TGAColor(255, 255, 255)*intensity; // well duh
		return false;                              // no, we do not discard this pixel
	}
};


struct GouraudShader : public IShader 
{
	Vec3f varying_intensity; // written by vertex shader, read by fragment shader

	virtual Vec4f vertex(int iface, int nthvert) {
		Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
		gl_Vertex = Viewport * Projection * View * World *gl_Vertex;     // transform it to screen coordinates
		varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
		return gl_Vertex;
	}

	virtual bool fragment(Vec3f bar, TGAColor &color) {
		float intensity = varying_intensity * bar;   // interpolate intensity for the current pixel
		color = TGAColor(255, 255, 255)*intensity; // well duh
		return false;                              // no, we do not discard this pixel
	}
};

int main()
{
	model = new Model("Models/african_head.obj");

	lookat(eye, center, up);
	viewport(width / 8, height / 8, width * 3 / 4, height * 3 / 4);
	projection(-1.f / (eye - center).norm());
	zaWarudo(yaw, pitch, roll, scalingFactors, position);
	light_dir.normalize();

	TGAImage image(width, height, TGAImage::RGB);
	TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);


	SpecularWithNormalMappingShader shader;
	shader.uniform_M = View*World;
	shader.uniform_MIT = (View * World).invert_transpose();
	for (int i = 0; i < model->nfaces(); i++) {
		Vec4f screen_coords[3];
		for (int j = 0; j < 3; j++) {
			screen_coords[j] = shader.vertex(i, j);
		}
		triangle(screen_coords, shader, image, zbuffer);
	}

	image.flip_vertically(); // to place the origin in the bottom left corner of the image
	zbuffer.flip_vertically();
	image.write_tga_file("output.tga");
	zbuffer.write_tga_file("zbuffer.tga");

	delete model;
	return 0;
}