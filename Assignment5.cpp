#define NUM_OPENGL_LIGHTS 8

#include <iostream>
#include <fstream>
#include <string>
#include <GL/glut.h>
#include <GL/glui.h>
#include "Shape.h"
#include "Cube.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Sphere.h"
#include "SceneParser.h"
#include "Camera.h"
#include "flatten.h"
#include "raycasting.h"
#include "ppm.h"

using namespace std;

int recursiveSteps = 1;
/** These are the live variables passed into GLUI ***/
int  isectOnly = 1;

int	 camRotU = 0;
int	 camRotV = 0;
int	 camRotW = 0;
int  viewAngle = 45;
float eyeX = 2;
float eyeY = 2;
float eyeZ = 2;
float lookX = -2;
float lookY = -2;
float lookZ = -2;

/** These are GLUI control panel objects ***/
int  main_window;
string filenamePath = "data/tests/earthcube.xml";
GLUI_EditText* filenameTextField = NULL;
GLubyte* pixels = NULL;
int pixelWidth = 0, pixelHeight = 0;
int screenWidth = 0, screenHeight = 0;

/** these are the global variables used for rendering **/
Cube* cube = new Cube();
Cylinder* cylinder = new Cylinder();
Cone* cone = new Cone();
Sphere* sphere = new Sphere();
SceneParser* parser = NULL;
Camera* camera = new Camera();

SceneLightData *light_data = NULL;
SceneGlobalData global_data;
int num_lights;

void setupCamera();
void updateCamera();
SceneColor getPixelColor(Point p, Vector ray, Vector lookV, int recursiveDepth, double curr_ks);

void setPixel(GLubyte* buf, int x, int y, int r, int g, int b) {
	buf[(y*pixelWidth + x) * 3 + 0] = (GLubyte)r;
	buf[(y*pixelWidth + x) * 3 + 1] = (GLubyte)g;
	buf[(y*pixelWidth + x) * 3 + 2] = (GLubyte)b;
}

void callback_start(int id) {
	cout << "start button clicked!" << endl;

	if (parser == NULL) {
		cout << "no scene loaded yet" << endl;
		return;
	}

	pixelWidth = screenWidth;
	pixelHeight = screenHeight;

	updateCamera();

	if (pixels != NULL) {
		delete pixels;
	}
	pixels = new GLubyte[pixelWidth  * pixelHeight * 3];
	memset(pixels, 0, pixelWidth  * pixelHeight * 3);

	cout << "(w, h): " << pixelWidth << ", " << pixelHeight << endl;
	cout << "Recursion: " << recursiveSteps << endl;
	
	Vector ray;

	// initialization of our stuff
	bool success;
    init_raycast(cube, cylinder, cone, sphere);
	success = parser->getGlobalData(global_data);

	num_lights = parser->getNumLights();
	if (light_data != NULL) {
		delete[] light_data;
	}
	light_data = new SceneLightData[num_lights];
	for (int i = 0; i < num_lights; ++i) {
		success = parser->getLightData(i, light_data[i]);
		if (!success) {
			std::cout << "could not load light on index = " << i << " with numLights = " << num_lights
				  << std::endl;
			exit(1);
		}
	}


	Point eyePoint = camera->GetEyePoint();
	Vector lookV = camera->GetLookVector();
	
	for (int c = 0; c < pixelWidth; c++) {
		for (int r = 0; r < pixelHeight; r++) {
			Vector ray = generateRay(eyePoint, c, r, pixelWidth, pixelHeight, camera);

			SceneColor color = getPixelColor(eyePoint, ray, lookV, 1, 1.0);

			setPixel(pixels, c, r, color.channels[0] * 255, color.channels[1] * 255, color.channels[2] * 255);
		}
	}
	glutPostRedisplay();
}

SceneColor getPixelColor(Point p, Vector ray, Vector lookV, int recursiveDepth, double curr_ks)
{
	SceneColor color;
	color.channels[0] = 0.0;
	color.channels[1] = 0.0;
	color.channels[2] = 0.0;

	if ((recursiveDepth > recursiveSteps) || (curr_ks < 0.001)) {
		return color;
	}

	double t = -1.0;
	double min_t = t;
	double shadow_t;
	bool shadow;
	FlatNode *iter = NULL;
	FlatInfo min_info;

	double kd = global_data.kd;
	double ks = global_data.ks;

	iter = flatRoot;
	t = -1.0;
	min_t = t;
	while (iter != NULL) {
		t = intersect(&(iter->info), p, ray);
		if (t > EPSILON && (t < min_t || min_t == -1.0)) {
			min_t = t;
			min_info = iter->info;
		}
		iter = iter -> next;
	}
	if (min_t > 0) {
		if (isectOnly) {
			color.channels[0] = 1.0;
			color.channels[1] = 1.0;
			color.channels[2] = 1.0;
		} else {
			SceneLightData light;
			Point intersect_p = p + min_t * ray;

			Matrix mt = min_info.transformation;
			Matrix mt_inv = invert(mt);

			Point p_obj = mt_inv * p;
			Vector ray_obj = mt_inv * ray;
			Point intersect_obj = mt_inv * intersect_p;

			Vector normal_obj = getNormal(p_obj, ray_obj, min_t, min_info);
			Vector normal_world = transpose(mt_inv) * normal_obj;
			normal_world.normalize();


			color = getAmbientContribution(global_data.ka, min_info);

			ScenePrimitive *prim = min_info.primitives[min_info.intersectIndex];
			SceneMaterial mat = prim->material;
			SceneFileMap *tex = mat.textureMap;

			SceneColor textureContrib;
			if (tex->isUsed) {
				TexInfo *tex_info = findTexture(tex->filename);
				if (tex_info == NULL) {
					std::cout << "tex_ppm pointer is NULL" << std::endl;
				}
				textureContrib = getTextureContrib(tex_info, intersect_obj, min_info);
			} else {
				textureContrib.channels[0] = 0.0;
				textureContrib.channels[1] = 0.0;
				textureContrib.channels[2] = 0.0;

			}
			for (int i = 0; i < num_lights; ++i) {
				light = light_data[i];

				Vector light_dir = calculateLightDir(&light, intersect_p);
				light_dir.normalize();
	
				// check all objects to see if blocking light
				shadow = false;
				iter = flatRoot;
				shadow_t = -1.0;
				Point new_intersect_p = intersect_p;
				new_intersect_p[0] += EPSILON;
				new_intersect_p[1] += EPSILON;
				new_intersect_p[2] += EPSILON;
				while (iter != NULL) {
					shadow_t = intersect(&(iter->info), new_intersect_p, -light_dir);
					if (shadow_t > EPSILON) {
						shadow = true;
						break;
					}
					iter = iter -> next;
				}

			    SceneColor diffuseContrib = getDiffuseContrib(kd, &(light), normal_world, light_dir,
			    											  min_info, intersect_obj, textureContrib, shadow);
			    SceneColor specularContrib = getSpecularContrib(ks, &(light), normal_world, light_dir,
			    										        ray, min_info, shadow);
                
            /*
                specularContrib.channels[0] = 0.0;
                specularContrib.channels[1] = 0.0;
                specularContrib.channels[2] = 0.0;
*/

			    color.channels[0] += (diffuseContrib.channels[0] + specularContrib.channels[0]);
			    color.channels[1] += (diffuseContrib.channels[1] + specularContrib.channels[1]);
			    color.channels[2] += (diffuseContrib.channels[2] + specularContrib.channels[2]);
			} 


			Vector reflected = getReflectedVector(ray, normal_world);
			reflected.normalize();

			SceneColor reflectiveContrib = 
				getPixelColor(intersect_p, reflected, lookV, recursiveDepth + 1, curr_ks * ks);
			reflectiveContrib.channels[0] *= (ks * prim->material.cReflective.channels[0]);
			reflectiveContrib.channels[1] *= (ks * prim->material.cReflective.channels[1]);
			reflectiveContrib.channels[2] *= (ks * prim->material.cReflective.channels[2]);	

			color.channels[0] += reflectiveContrib.channels[0];
			color.channels[1] += reflectiveContrib.channels[1];
			color.channels[2] += reflectiveContrib.channels[2];


			if (color.channels[0] > 1.0)
				color.channels[0] = 1.0;
			if (color.channels[1] > 1.0)
				color.channels[1] = 1.0;
			if (color.channels[2] > 1.0)
				color.channels[2] = 1.0;
		}
	} else {
		color.channels[0] = 0.0;
		color.channels[1] = 0.0;
		color.channels[2] = 0.0;
	}

	return color;
}



void callback_load(int id) {
	char curDirName [2048];
	if (filenameTextField == NULL) {
		return;
	}
	printf ("%s\n", filenameTextField->get_text());

	if (parser != NULL) {
		delete parser;
	}
	parser = new SceneParser (filenamePath);
	cout << "success? " << parser->parse() << endl;

	//FLATTEN HERE
	SceneNode *root = parser->getRootNode();
	if (flatRoot != NULL) {
		deleteList();
	}
	flatten(root, true);
	if (texRoot != NULL) {
		deleteTextureList();
	}
	flattenTextures();

	setupCamera();
}


/***************************************** myGlutIdle() ***********/

void myGlutIdle(void)
{
	/* According to the GLUT specification, the current window is
	undefined during an idle callback.  So we need to explicitly change
	it if necessary */
	if (glutGetWindow() != main_window)
		glutSetWindow(main_window);

	glutPostRedisplay();
}


/**************************************** myGlutReshape() *************/

void myGlutReshape(int x, int y)
{
	float xy_aspect;

	xy_aspect = (float)x / (float)y;
	glViewport(0, 0, x, y);
	camera->SetScreenSize(x, y);

	screenWidth = x;
	screenHeight = y;

	glutPostRedisplay();
}


/***************************************** setupCamera() *****************/
void setupCamera()
{
	SceneCameraData cameraData;
	parser->getCameraData(cameraData);

	camera->Reset();
	camera->SetViewAngle(cameraData.heightAngle);
	if (cameraData.isDir == true) {
		camera->Orient(cameraData.pos, cameraData.look, cameraData.up);
	}
	else {
        cout << "cameraData.pos = " << cameraData.pos[0] << " " << cameraData.pos[1] << " " << cameraData.pos[2] << endl;
        cout << "cameraData.lookAt = " << cameraData.lookAt[0] << " " << cameraData.lookAt[1] << " " << cameraData.lookAt[2] << endl;
		camera->Orient(cameraData.pos, cameraData.lookAt, cameraData.up);
	}

	viewAngle = camera->GetViewAngle();
	Point eyeP = camera->GetEyePoint();
	Vector lookV = camera->GetLookVector();
    lookV.normalize();
    cout << "lookv = " << lookV[0] << " " << lookV[1] << " " << lookV[2] << endl;
    cout << "eye point = " << eyeP[0] << " " << eyeP[1] << " " << eyeP[2] << endl;
	eyeX = eyeP[0];
	eyeY = eyeP[1];
	eyeZ = eyeP[2];
	lookX = lookV[0];
	lookY = lookV[1];
	lookZ = lookV[2];
	camRotU = 0;
	camRotV = 0;
	camRotW = 0;
	GLUI_Master.sync_live_all();
}

void updateCamera()
{
	camera->Reset();

	Point guiEye (eyeX, eyeY, eyeZ);
	Point guiLook(lookX, lookY, lookZ);
	camera->SetViewAngle(viewAngle);
	Vector upV = camera->GetUpVector();
	//camera->Orient(guiEye, guiLook, camera->GetUpVector());
	camera->Orient(guiEye, guiLook, upV);
	camera->RotateU(camRotU);
	camera->RotateV(camRotV);
	camera->RotateW(camRotW);
}

/***************************************** myGlutDisplay() *****************/

void myGlutDisplay(void)
{
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT);

	if (parser == NULL) {
		return;
	}

	if (pixels == NULL) {
		return;
	}
	
	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
	glDrawPixels(pixelWidth, pixelHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels);
	glutSwapBuffers();
}

void onExit()
{
	delete cube;
	delete cylinder;
	delete cone;
	delete sphere;
	delete camera;
	if (parser != NULL) {
		delete parser;
	}
	if (pixels != NULL) {
		delete pixels;
	}
}

/**************************************** main() ********************/

int main(int argc, char* argv[])
{
        fprintf(stderr, "Hello world!\n");
	atexit(onExit);

	/****************************************/
	/*   Initialize GLUT and create window  */
	/****************************************/

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(500, 500);

	main_window = glutCreateWindow("COMP 175 Assignment 4");
	glutDisplayFunc(myGlutDisplay);
	glutReshapeFunc(myGlutReshape);

	/****************************************/
	/*         Here's the GLUI code         */
	/****************************************/

	GLUI* glui = GLUI_Master.create_glui("GLUI");

	filenameTextField = new GLUI_EditText( glui, "Filename:", filenamePath);
	filenameTextField->set_w(300);
	glui->add_button("Load", 0, callback_load);
	glui->add_button("Start!", 0, callback_start);
	glui->add_checkbox("Isect Only", &isectOnly);

	GLUI_Panel *recursive_panel = glui->add_panel("Recursion");
	(new GLUI_Spinner(recursive_panel, "Recursive Depth:", &recursiveSteps))
		->set_int_limits(1, 25);
	glui->add_column_to_panel(recursive_panel, true);
	
	GLUI_Panel *camera_panel = glui->add_panel("Camera");
	(new GLUI_Spinner(camera_panel, "RotateV:", &camRotV))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "RotateU:", &camRotU))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "RotateW:", &camRotW))
		->set_int_limits(-179, 179);
	(new GLUI_Spinner(camera_panel, "Angle:", &viewAngle))
		->set_int_limits(1, 179);

	glui->add_column_to_panel(camera_panel, true);

	GLUI_Spinner* eyex_widget = glui->add_spinner_to_panel(camera_panel, "EyeX:", GLUI_SPINNER_FLOAT, &eyeX);
	eyex_widget->set_float_limits(-10, 10);
	GLUI_Spinner* eyey_widget = glui->add_spinner_to_panel(camera_panel, "EyeY:", GLUI_SPINNER_FLOAT, &eyeY);
	eyey_widget->set_float_limits(-10, 10);
	GLUI_Spinner* eyez_widget = glui->add_spinner_to_panel(camera_panel, "EyeZ:", GLUI_SPINNER_FLOAT, &eyeZ);
	eyez_widget->set_float_limits(-10, 10);

	GLUI_Spinner* lookx_widget = glui->add_spinner_to_panel(camera_panel, "LookX:", GLUI_SPINNER_FLOAT, &lookX);
	lookx_widget->set_float_limits(-10, 10);
	GLUI_Spinner* looky_widget = glui->add_spinner_to_panel(camera_panel, "LookY:", GLUI_SPINNER_FLOAT, &lookY);
	looky_widget->set_float_limits(-10, 10);
	GLUI_Spinner* lookz_widget = glui->add_spinner_to_panel(camera_panel, "LookZ:", GLUI_SPINNER_FLOAT, &lookZ);
	lookz_widget->set_float_limits(-10, 10);

	glui->add_button("Quit", 0, (GLUI_Update_CB)exit);

	glui->set_main_gfx_window(main_window);

	/* We register the idle callback with GLUI, *not* with GLUT */
	GLUI_Master.set_glutIdleFunc(myGlutIdle);

	glutMainLoop();

	return EXIT_SUCCESS;
}



