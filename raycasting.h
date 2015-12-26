// Jake Mingolla & Casey Gowrie
// April 2015
// Comp 175 Graphics
// Remco

#ifndef RAYCASTING_H
#define RAYCASTING_H

#include <cmath>
#include "flatten.h"
#include "Camera.h"
#include "Shape.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "SceneData.h"
#include "ppm.h"

Cube *R_cube = NULL;
Cylinder *R_cylinder = NULL;
Cone *R_cone = NULL;
Sphere *R_sphere = NULL;

void init_raycast(Cube *cube, Cylinder *cylinder, Cone *cone, Sphere *sphere) {
    R_cube = cube;
    R_cylinder = cylinder;
    R_cone = cone;
    R_sphere = sphere;
}


Vector generateRay(Point P, int c, int r, int width, int height, Camera *camera) {
	Vector ray;
	double nearPlane, theta, aspectRatio, near;
	double fp_height, fp_width;
	double a, b;
	Point Q, S;
	Vector lookV, upV, u, v, w;

	lookV = camera->GetLookVector();
	lookV.normalize();
	upV = camera->GetUpVector();
	upV.normalize();
	near = camera->GetNearPlane();
	Q = P + (lookV * near);

	w = -1 * lookV;
	w = (1.0 / w.length()) * w;
	u = cross(upV, w);
	u = u * (1.0 / u.length());
	v = cross(w, u);

	theta = camera->GetViewAngle();
	aspectRatio = camera->GetScreenWidthRatio();
	fp_height = near * tan(DEG_TO_RAD(theta / 2));
	fp_width = fp_height * aspectRatio;

	a = (-1 * fp_width) + ((2 * fp_width) * ((float)c / (float)width));
	b = (-1 * fp_height) + ((2 * fp_height) * ((float)r / (float)height));

	S = Q;
	S = S + (a * u);
	S = S + (b * v);

	ray = S - P;
	ray.normalize();

	return ray;
}

double intersect(FlatInfo *info, Point p, Vector ray) {
	Matrix inv_trans = invert(info->transformation);
	Point eye_obj = inv_trans * p;
	Vector ray_obj = inv_trans * ray;

	double t = -1.0;
	double min_t = -1.0;
	
	int len = info->primitives.size();
	ScenePrimitive *prim;
	
	for (int i = 0; i < len; ++i) {
		prim = info->primitives[i];
		switch(prim->type) {
			case SHAPE_CUBE:
				// CUBE
				t = R_cube->intersect(eye_obj, ray_obj);
				break;
			case SHAPE_CYLINDER:
				//cylinder
				t = R_cylinder->intersect(eye_obj, ray_obj);
				break;
			case SHAPE_CONE:
				//cone
				t = R_cone->intersect(eye_obj, ray_obj);
				break;
			case SHAPE_SPHERE:
				//sphere
				t = R_sphere->intersect(eye_obj, ray_obj);
				break;
		}
		if (t > 0 && (t < min_t || min_t == -1.0)) {
			min_t = t;
			info->intersectIndex = i;
		}
	}
	
	return min_t;
}

Vector getNormal(Point p, Vector ray, double t, FlatInfo info) {
	ScenePrimitive *prim = info.primitives[info.intersectIndex];
	Vector n;	
	switch(prim->type) {
			case SHAPE_CUBE:
				// CUBE
				n = R_cube->getIsectNormal(p, ray, t);
				break;
			case SHAPE_CYLINDER:
				//cylinder
				n = R_cylinder->getIsectNormal(p, ray, t);
				break;
			case SHAPE_CONE:
				//cone
				n = R_cone->getIsectNormal(p, ray, t);
				break;
			case SHAPE_SPHERE:
				//sphere
				n = R_sphere->getIsectNormal(p, ray, t);
				break;
		}
	return n;
}

SceneColor getAmbientContribution(float ka, FlatInfo info) {
	ScenePrimitive *prim = info.primitives[info.intersectIndex];	
	SceneMaterial mat = prim->material;
	SceneColor ambient, cA = mat.cAmbient;

	ambient.r = ka * cA.channels[0];
	ambient.g = ka * cA.channels[1];
	ambient.b = ka * cA.channels[2];
	return ambient;
}


Point mapToUnitSquare(ScenePrimitive *prim, Point p)
{
	Point unit_square_p = Point(0, 0, 0);

	switch(prim->type) {
			case SHAPE_CUBE:
				// CUBE
				unit_square_p = R_cube->mapToSquare(p);
				break;
			case SHAPE_CYLINDER:
				unit_square_p = R_cylinder->mapToSquare(p);
				//cylinder
				break;
			case SHAPE_CONE:
				//cone
				unit_square_p = R_cone->mapToSquare(p);
				break;
			case SHAPE_SPHERE:
				//sphere
				unit_square_p = R_sphere->mapToSquare(p);
				break;
	}

	return unit_square_p;
}

Point unitSquareToTexture(Point us, int w, int h, float u, float v) {
	Point tex_p;

	tex_p[0] = (int)(us[0] * (float) w * u) % w;
	tex_p[1] = (int)(us[1] * (float) h * v) % h;

	return tex_p;
}

SceneColor getTextureContrib(TexInfo *tex_info, Point intersect_obj, FlatInfo info) {

	Point unit, tex_p; //only 2d point
    int tex_width, tex_height, tex_len;

	SceneColor tex_color;
	tex_color.channels[0] = 0.0;
	tex_color.channels[1] = 0.0;
	tex_color.channels[2] = 0.0;
	ScenePrimitive *prim = info.primitives[info.intersectIndex];
	SceneMaterial mat = prim->material;
	SceneFileMap *tex = mat.textureMap;
	if (tex->isUsed) {
		unit = mapToUnitSquare(prim, intersect_obj);

		tex_width = tex_info->width;
		tex_height = tex_info->height;
        tex_len = tex_info->len;
		tex_p = unitSquareToTexture(unit, tex_width, tex_height, tex->repeatU, tex->repeatV);
		char *tex_array = tex_info->arr;

		int tex_x = (int) tex_p[0];
		int tex_y = (int) tex_p[1];
		int tex_index_r = 3 * (tex_y * tex_width + tex_x);

		unsigned int r = tex_array[tex_index_r + 0] + 0;
		unsigned int g = tex_array[tex_index_r + 1] + 0;
		unsigned int b = tex_array[tex_index_r + 2] + 0;

        r = r << 24;
        r = r >> 24;
        g = g << 24;
        g = g >> 24;
        b = b << 24;
        b = b >> 24;


		tex_color.channels[0] = (float) r / 255.0;
		tex_color.channels[1] = (float) g / 255.0;
		tex_color.channels[2] = (float) b / 255.0;
	}
	return tex_color;
}

SceneColor getDiffuseContrib(float kd, SceneLightData *light, Vector normal, Vector light_dir, 
							FlatInfo info, Point intersect_obj, SceneColor texture, bool shadow) 
{
	ScenePrimitive *prim = info.primitives[info.intersectIndex];
	SceneMaterial mat = prim->material;
	SceneColor diffuse, cD = mat.cDiffuse, lC = light->color;
	double dot_prod = dot(-light_dir, normal);
	double tex_weight = mat.blend;
	double dif_weight = 1.0 - tex_weight;

	if (dot_prod > 0.0 && !shadow) {
        // blend after diffuse
		diffuse.channels[0] = kd * cD.channels[0] *
							  lC.channels[0] * dot_prod;
		diffuse.channels[1] = kd * cD.channels[1] *
							  lC.channels[1] * dot_prod;
		diffuse.channels[2] = kd * cD.channels[2] *
							  lC.channels[2] * dot_prod;
	} else {
		diffuse.channels[0] = 0.0;
		diffuse.channels[1] = 0.0;
		diffuse.channels[2] = 0.0;
	}
	//cout << "kd = " << kd << endl;
	if (texture.channels[0] > 1.0) {
		texture.channels[0] = 1.0;
	}
	if (texture.channels[1] > 1.0) {
		texture.channels[1] = 1.0;
	}
	if (texture.channels[2] > 1.0) {
		texture.channels[2] = 1.0;
	}
	diffuse.channels[0] = diffuse.channels[0] * dif_weight + (texture.channels[0] * tex_weight);
    diffuse.channels[1] = diffuse.channels[1] * dif_weight + (texture.channels[1] * tex_weight);
    diffuse.channels[2] = diffuse.channels[2] * dif_weight + (texture.channels[2] * tex_weight);

	return diffuse;
}

Vector getReflectedVector(Vector v, Vector n) {
	return v - (2 * dot(v, n) * n);
}

SceneColor getSpecularContrib(float ks, SceneLightData *light, Vector normal, Vector light_dir,
							  Vector lookV, FlatInfo info, bool shadow) {
	ScenePrimitive *prim = info.primitives[info.intersectIndex];
	SceneMaterial mat = prim->material;
	SceneColor specular, cS = mat.cSpecular, lC = light->color;

	//cout << "light_dir: " << light_dir[0] << " " << light_dir[1] << " " << light_dir[2] << endl;

	Vector reflected = getReflectedVector(light_dir, normal);
	reflected.normalize();
	lookV.normalize();

	double dot_product = dot(-lookV, reflected);
	double dot_prod_to_f = pow(dot_product, mat.shininess);
	// 2 above

	//if (dot_product > 0.0 && !shadow) {
	// 1 above
	//cout << "SHINE: " << mat.shininess << endl;
	if (dot_product > 0.0 && !shadow) {
		specular.channels[0] = ks * cS.channels[0] * dot_prod_to_f * lC.channels[0];
		specular.channels[1] = ks * cS.channels[1] * dot_prod_to_f * lC.channels[1];
		specular.channels[2] = ks * cS.channels[2] * dot_prod_to_f * lC.channels[2];
	} else if (IN_RANGE(dot_prod_to_f, 1.0) && !shadow) {
		specular.channels[0] = ks * cS.channels[0] * lC.channels[0];
		specular.channels[1] = ks * cS.channels[1] * lC.channels[1];
		specular.channels[2] = ks * cS.channels[2] * lC.channels[2];
	} else {
		specular.channels[0] = 0.0;
		specular.channels[1] = 0.0;
		specular.channels[2] = 0.0;
	}
	return specular;
}

Vector calculateLightDir(SceneLightData *light, Point intersect) {
	Vector dir = Vector(0, 0, 0);

	switch (light->type) {
		case LIGHT_POINT:
			dir = intersect - light->pos;
			break;
		case LIGHT_DIRECTIONAL:
			dir = light->dir;
			break;
		case LIGHT_SPOT:
		// how? need distance to lights??
			break;
		case LIGHT_AREA:
		// how? need locations from light??
			break;
	}

	return dir;
}
	

#endif


