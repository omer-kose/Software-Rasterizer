#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include "GeneralGraphicsUtilites.h"

const double pi = 3.1415926535897932385;


Matrix View;
Matrix Viewport;
Matrix Projection;
Matrix World;

IShader::~IShader() {}


inline double degrees_to_radians(double degrees)
{
	return degrees * pi / 180.0;
}

/*XYZ respectively*/
Matrix scale(Vec3f scalingFactors)
{
	Matrix scaling = Matrix::identity();
	for (int i = 0; i < 3; ++i)
	{
		scaling[i][i] = scalingFactors[i];
	}

	return scaling;
}
Matrix translate(Vec3f position)
{
	Matrix translation = Matrix::identity();
	for (int i = 0; i < 3; ++i)
	{
		translation[i][3] = position[i];
	}
	return translation;
}
Matrix rotationX(double angle)
{
	Matrix rotation = Matrix::identity();
	angle = degrees_to_radians(angle);
	double sine = sin(angle);
	double cosine = cos(angle);
	rotation[1][1] = cosine;
	rotation[2][1] = sine;
	rotation[1][2] = -sine;
	rotation[2][2] = cosine;
	
	return rotation;
}

Matrix rotationY(double angle)
{
	Matrix rotation = Matrix::identity();
	angle = degrees_to_radians(angle);
	double sine = sin(angle);
	double cosine = cos(angle);
	rotation[0][0] = cosine;
	rotation[2][0] = -sine;
	rotation[0][2] = sine;
	rotation[2][2] = cosine;

	return rotation;
}

Matrix rotationZ(double angle)
{
	Matrix rotation = Matrix::identity();
	angle = degrees_to_radians(angle);
	double sine = sin(angle);
	double cosine = cos(angle);
	rotation[0][0] = cosine;
	rotation[1][0] = sine;
	rotation[0][1] = -sine;
	rotation[1][1] = cosine;

	return rotation;
}


/*Note that this rasterizer can only render one model at each run*/
/*It sets the world matrix of our model*/
/*Yeah this is a JOJO reference*/
void zaWarudo(double yaw, double pitch, double roll, Vec3f scalingFactors, Vec3f position)
{
	World = Matrix::identity();
	/*Rotation*/
	Matrix mRotationX, mRotationY, mRotationZ, mRotation;
	mRotationX = rotationX(pitch);
	mRotationY = rotationY(yaw);
	mRotationZ = rotationZ(roll);
	mRotation = mRotationX * mRotationY * mRotationZ;

	/*Scaling and Translation*/
	Matrix scaling = scale(scalingFactors);
	Matrix translation = translate(position);

	World = translation * scaling * mRotation * World;
}

void viewport(int x, int y, int w, int h)
{
	Viewport = Matrix::identity();
	Viewport[0][3] = x + w / 2.f;
	Viewport[1][3] = y + h / 2.f;
	Viewport[2][3] = 255.f / 2.f;
	Viewport[0][0] = w / 2.f;
	Viewport[1][1] = h / 2.f;
	Viewport[2][2] = 255.f / 2.f;
}

void projection(float coeff) 
{
	Projection = Matrix::identity();
	Projection[3][2] = coeff;
}

void lookat(Vec3f eye, Vec3f center, Vec3f up)
{
	Vec3f z = (eye - center).normalize();
	Vec3f x = cross(up, z).normalize();
	Vec3f y = cross(z, x).normalize();
	View = Matrix::identity();
	for (int i = 0; i < 3; i++) {
		View[0][i] = x[i];
		View[1][i] = y[i];
		View[2][i] = z[i];
		View[i][3] = -center[i];
	}
}

Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P)
{
	Vec3f s[2];
	for (int i = 2; i--; ) {
		s[i][0] = C[i] - A[i];
		s[i][1] = B[i] - A[i];
		s[i][2] = A[i] - P[i];
	}
	Vec3f u = cross(s[0], s[1]);
	if (std::abs(u[2]) > 1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
		return Vec3f(1.f - (u.x + u.y) / u.z, u.y / u.z, u.x / u.z);
	return Vec3f(-1, 1, 1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

/*For now, I will stick to this triangle function*/
void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer) {
	Vec2f bboxmin(std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {
			bboxmin[j] = std::min(bboxmin[j], pts[i][j] / pts[i][3]);
			bboxmax[j] = std::max(bboxmax[j], pts[i][j] / pts[i][3]);
		}
	}
	Vec2i P;
	TGAColor color;
	for (P.x = bboxmin.x; P.x <= bboxmax.x; P.x++) {
		for (P.y = bboxmin.y; P.y <= bboxmax.y; P.y++) {
			Vec3f c = barycentric(proj<2>(pts[0] / pts[0][3]), proj<2>(pts[1] / pts[1][3]), proj<2>(pts[2] / pts[2][3]), proj<2>(P));
			float z = pts[0][2] * c.x + pts[1][2] * c.y + pts[2][2] * c.z;
			float w = pts[0][3] * c.x + pts[1][3] * c.y + pts[2][3] * c.z;
			int frag_depth = std::max(0, std::min(255, int(z / w + .5)));
			if (c.x < 0 || c.y < 0 || c.z<0 || zbuffer.get(P.x, P.y)[0]>frag_depth) continue;
			bool discard = shader.fragment(c, color);
			if (!discard) {
				zbuffer.set(P.x, P.y, TGAColor(frag_depth));
				image.set(P.x, P.y, color);
			}
		}
	}
}

/*Will be worked out later, barycentric approach will be aplied*/
//void triangleMine(Vec4f *vertices, IShader &shader, TGAImage &image, TGAImage &zbuffer)
//{
//	/*sort the vertices so that from up to bottom the order is v0 v1 v2*/
//	if (vertices[0][1] > vertices[1][1]) std::swap(vertices[0], vertices[1]);
//	if (vertices[0][1] > vertices[2][1]) std::swap(vertices[0], vertices[2]);
//	if (vertices[1][1] > vertices[2][1]) std::swap(vertices[1], vertices[2]);
//
//
//	/*Find VPrime that is the vertex that lies in the same y with middle vertex(v1)*/
//	/*Find by interpolating from v0 to v2 using alpha as v1's y*/
//	/*Using floats to avoid preicison errors, while putting pixels points will be casted to ints*/
//	float yPrime = vertices[1][1];
//	float alphaPrime = (yPrime - vertices[2][1]) / (float)(vertices[0][1] - vertices[2][1]);
//	float xPrime = vertices[2][0] + alphaPrime * (vertices[0][0] - vertices[2][0]);
//
//	/*iterate from y0 to y1 and fill upper half*/
//	for (int y = vertices[2][1]; y > vertices[1][1]; y--)
//	{
//		float alpha = (y - vertices[2][1]) / (float)(vertices[1][1] - vertices[2][1]);
//		/*Determine left boundary*/
//		int xLeft = vertices[2][0] + alpha * (vertices[1][0] - vertices[2][0]);
//		/*Determine right boundary*/
//		int xRight = vertices[2][0] + alpha * (xPrime - vertices[2][0]);
//
//		
//		if (xLeft > xRight)
//		{
//			std::swap(xLeft, xRight);
//		}
//
//		/*Fill the pixels*/
//		for (int x = xLeft; x <= xRight; x++)
//		{
//			Vec3f P(x, y, 0);
//			Vec3f bc_screen = barycentric(vertices[0], vertices[1], vertices[2], P);
//			for (int i = 0; i < 3; i++)  P.z += vertices[i][2] * bc_screen[i];
//			if (zBuffer[int(P[1] * width + P[0])] < P.z)
//			{
//				zBuffer[int(P[1]*width + P[0])] = P.z;
//				image.set(x, y, color);
//			}
//			
//		}
//
//	}
//
//	/*Draw the lower half*/
//	for (int y = vertices[1][1]; y > vertices[0][1]; y--)
//	{
//		float alpha = (y - vertices[1][1]) / (float)(vertices[0][1] - vertices[1][1]);
//		/*Determine left boundary*/
//		int xLeft = vertices[1][0] + alpha * (vertices[0][0] - vertices[1][0]);
//		/*Determine right boundary*/
//		int xRight = xPrime + alpha * (vertices[0][0] - xPrime);
//
//
//		if (xLeft > xRight)
//		{
//			std::swap(xLeft, xRight);
//		}
//
//		/*Fill the pixels*/
//		for (int x = xLeft; x <= xRight; x++)
//		{
//			Vec3f P(x, y, 0);
//			Vec3f bc_screen = barycentric(vertices[0], vertices[1], vertices[2], P);
//			for (int i = 0; i < 3; i++)  P.z += vertices[i][2] * bc_screen[i];
//			if (zBuffer[int(P[1] * width + P[0])] < P.z)
//			{
//				zBuffer[int(P[1]*width + P[0])] = P.z;
//				image.set(x, y, color);
//			}
//		}
//
//	}
//}
