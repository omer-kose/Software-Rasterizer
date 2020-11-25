#pragma once
#include "tgaimage.h"
#include "geometry.h"

extern Matrix View;
extern Matrix Viewport;
extern Matrix Projection;
extern Matrix World;

const float depth = 2000.f;


void viewport(int x, int y, int w, int h);
void projection(float coeff = 0.f); // coeff = -1/c
void lookat(Vec3f eye, Vec3f center, Vec3f up);
void zaWarudo(double yaw, double pitch, double roll, Vec3f scalingFactors, Vec3f position);

struct IShader {
	virtual ~IShader();
	virtual Vec4f vertex(int iface, int nthvert) = 0;
	virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};

void triangle(Vec4f *vertices, IShader &shader, TGAImage &image, TGAImage &zbuffer);
