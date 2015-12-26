/*
Casey Gowrie and Jake Mingolla
March, 2015
Comp 175 -- Computer Graphics
*/

#ifndef CAMERA_H
#define CAMERA_H

#include "Algebra.h"

class Camera {
public:
	Camera();
	~Camera();

	void Reset();

	// in every one of these set bool to true if different
	void Orient(Point& new_eye, Point& new_focus, Vector& new_up);
	void Orient(Point& new_eye, Vector& new_look, Vector& new_up);
	void SetViewAngle (double new_viewAngle);
	void SetNearPlane (double new_nearPlane);
	void SetFarPlane (double new_farPlane);
	void SetScreenSize (int new_screenWidth, int new_screenHeight);

	// if changed recalculate, set changed to false
	Matrix GetProjectionMatrix();
	Matrix GetModelViewMatrix();

	void SetAngleU(double angle);
	void SetAngleV(double angle);
	void SetAngleW(double angle);

	void RotateV(double angle);
	void RotateU(double angle);
	void RotateW(double angle);
	void Rotate(Point p, Vector axis, double degree);
	void Translate(const Vector &v);

	Point GetEyePoint();
	Vector GetLookVector();
	Vector GetUpVector();
	double GetViewAngle();
	double GetNearPlane();
	double GetFarPlane();
	int GetScreenWidth();
	int GetScreenHeight();

	double GetFilmPlanDepth();
	double GetScreenWidthRatio();

private:
	Point eye;
	Vector look, up;
	Vector u, v, w;
	double viewAngle, nearPlane, farPlane;
	int screenWidth, screenHeight;
	double angle_U, angle_V, angle_W;
	Matrix M_Translate, M_Rotate, M_Scale, M_Unhinge;
	Matrix Projection_Matrix, ModelView_Matrix;

	void updateModelViewMatrix();
	void updateProjectionMatrix();
	double calculateWidthAngle();

	void printMatrix(Matrix m);
	void printVector(Vector v);

	bool changed;

};
#endif


