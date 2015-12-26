/*
Casey Gowrie and Jake Mingolla
March, 2015
Comp 175 -- Computer Graphics
*/
#include "Camera.h"
#include <iostream>

Camera::Camera() {
	M_Translate = Matrix();
	M_Rotate = Matrix();
	M_Scale = Matrix();
	M_Unhinge = Matrix();
	ModelView_Matrix = Matrix();
	Projection_Matrix = Matrix();
	u = Vector();
	v = Vector();
	w = Vector();
}

Camera::~Camera() {
}

void Camera::Reset() {
	nearPlane = 0.01;
	farPlane = 30;
	viewAngle = 45;
}

void Camera::Orient(Point& new_eye, Point& new_focus, Vector& new_up) {
	// TODO, make look from point and eye
	eye = new_eye;
	look = new_focus - new_eye;
	up = new_up;
	w = -1 * look;
	w = (1.0 / w.length()) * w;
	u = cross(up, w);
	u = u * (1.0 / u.length());
	v = cross(w, u);
}



void Camera::Orient(Point& new_eye, Vector& new_look, Vector& new_up) {
		eye = new_eye;
		look = new_look;
		up = new_up;

		w = -1 * look;
		w = (1.0 / w.length()) * w;

		u = cross(up, w);
		u = u * (1.0 / u.length());

		v = cross(w, u);
}

double Camera::calculateWidthAngle() {
	return viewAngle * ((float) screenWidth / (float) screenHeight);
}

void Camera::updateProjectionMatrix() {
	double widthAngle = calculateWidthAngle();
	double widthRad = widthAngle * (PI / 180);
	double heightRad = viewAngle * (PI / 180);
	M_Scale = Matrix(1 / (tan(widthRad/2)*farPlane), 0, 						  0, 	   0,
					 0,							1/(tan(heightRad/2)*farPlane), 0, 	   0,
					 0, 						0, 						  1 / farPlane, 0,
					 0,							0, 						  0,       1);

	double c = -1 * nearPlane / farPlane;
	M_Unhinge = Matrix(1, 0, 0, 	     0,
					   0, 1, 0, 	     0,
					   0, 0, -1/(c+1),   c/(c+1),
					   0, 0, -1, 0);


	Projection_Matrix = (M_Unhinge) * (M_Scale);

}

Matrix Camera::GetProjectionMatrix() {
	updateProjectionMatrix();

	return Projection_Matrix;
}


void Camera::SetViewAngle (double new_viewAngle) {
	if (!IN_RANGE(new_viewAngle, viewAngle)) {
		viewAngle = new_viewAngle;
		changed = true;
	}
}

void Camera::SetNearPlane (double new_nearPlane) {
	if (!IN_RANGE(new_nearPlane, nearPlane)) {
		nearPlane = new_nearPlane;
		changed = true;
	}
}

void Camera::SetFarPlane (double new_farPlane) {
	if (!IN_RANGE(new_farPlane, farPlane)) {
		farPlane = new_farPlane;
		changed = true;
	}
}

void Camera::SetScreenSize (int new_screenWidth, int new_screenHeight) {
	if ((new_screenWidth != screenWidth) || (new_screenHeight != screenHeight)) {
		screenWidth = new_screenWidth;
		screenHeight = new_screenHeight;
		changed = true;
	}
}

void Camera::updateModelViewMatrix() {
	// Only changing three values, do we need to change to identity?
	Point p_n = Point();
	p_n[0] = eye[0] + (nearPlane * look[0]);
	p_n[1] = eye[1] + (nearPlane * look[1]);
	p_n[2] = eye[2] + (nearPlane * look[2]);


	M_Translate = Matrix(1, 		   0, 			0, 			-1 * eye[0],
						 0, 		   1, 			0, 			-1 * eye[1],
						 0, 		   0, 		    1, 			-1 * eye[2],
						 0, 0, 0, 1);



	M_Rotate = Matrix(u[0], u[1], u[2], 0,
				  	  v[0], v[1], v[2], 0,
				  	  w[0], w[1], w[2], 0,
				  	  0,   0,   0,   1);

	ModelView_Matrix = (M_Rotate) * (M_Translate);
}

// if changed recalculate, set changed to false
Matrix Camera::GetModelViewMatrix() {
	updateModelViewMatrix();

	return ModelView_Matrix;
}


void Camera::RotateU(double angle) {
	//SetAngleU(angle_U + angle);
	angle *= (PI / 180);
	Matrix M_R = rot_mat(u, angle);

	v = M_R * v;
	w = M_R * w;
	look = M_R * look;
}

// calculate new look and/or up
void Camera::SetAngleU(double angle) {
	// if (!IN_RANGE(angle_U, angle)) {
	// 	double rotation = angle - angle_U;

	// 	angle_U = angle;
	// 	changed = true;
	// }
	return;

}

void Camera::RotateV(double angle) {
	//SetAngleV(angle_V + angle);
	angle *= (PI / 180);
	Matrix M_R = rot_mat(v, angle);

	u = M_R * u;
	w = M_R * w;
	look = M_R * look;
}

void Camera::SetAngleV(double angle) {
	return;
	// if (!IN_RANGE(angle, angle_V)) {
	// 	angle_V = angle;
	// 	changed = true;
	// 	std::cout << "changed angle v to " << angle << std::endl;
	// }
}

void Camera::RotateW(double angle) {
	angle *= (PI / 180);

	Matrix M_R = rot_mat(w, angle);
	u = M_R * u;
	v = M_R * v;
	up = M_R * up;

}

void Camera::SetAngleW(double angle) {
	// if (!IN_RANGE(angle_W, angle)) {
	// 	angle_W = angle;
	// 	changed = true;
	// }
	return;
}

void Camera::Translate(const Vector &v) {
	eye = eye + v;
}


void Camera::Rotate(Point p, Vector axis, double degrees) {
	// use rot_mat func again?? then rotate all?
	Matrix M_R = rot_mat(p, axis, degrees);
	up = M_R * up;
	look = M_R * look;
	u = M_R * u;
	v = M_R * v;
	w = M_R * w;
}


Point Camera::GetEyePoint() {
	return eye;
}

Vector Camera::GetLookVector() {
	return look;
}

Vector Camera::GetUpVector() {
	return up;
}

double Camera::GetViewAngle() {
	return viewAngle;
}

double Camera::GetNearPlane() {
	return nearPlane;
}

double Camera::GetFarPlane() {
	return farPlane;
}

int Camera::GetScreenWidth() {
	return screenWidth;
}

int Camera::GetScreenHeight() {
	return screenHeight;
}


// ASK ABOUT THE BOTTOM TWO FUNCTIONS

double Camera::GetFilmPlanDepth() {
	return farPlane - nearPlane;
}

double Camera::GetScreenWidthRatio() {
	return screenWidth / (double) screenHeight;
}

void Camera::printMatrix(Matrix m) {
	m = transpose(m);
	std::cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << std::endl;
	std::cout << m[4] << " " << m[5] << " " << m[6] << " " << m[7] << std::endl;
	std::cout << m[8] << " " << m[9] << " " << m[10] << " " << m[11] << std::endl;
	std::cout << m[12] << " " << m[13] << " " << m[14] << " " << m[15] << std::endl;
}

void Camera::printVector(Vector v) {
	std::cout << v[0] << " " << v[1] << " " << v[2] << std::endl;
}

