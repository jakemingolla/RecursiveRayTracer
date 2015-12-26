/* Jake Mingolla & Casey Gowrie
 * Shape parent class for assignment 1
 * Comp 175 -- Remco Chang */

#ifndef SHAPE_H
#define SHAPE_H

#include <GL/glui.h>
#include "Algebra.h"
#include "VertexList.h"
#include <iostream>

enum IntersectPlane {
    X = 0,
    Y = 1,
    Z = 2,
};


class Shape {
public:
	Shape() {
		x_dist = y_dist = z_dist = 1.0;
        m_segmentsX = m_segmentsY = 0;
		crossSections = NULL;
        numCrossSections = 0;
	};
    /* Since crossSections is the only dynamically allocated part of a 
 *     shape child class that is shared among all children, it is 
 *     deleted here. */
	~Shape() {
		if (crossSections != NULL) {
			delete[] crossSections;
		}
	};
    /* Given an x and y values for number of segments, if there is a
 *     change in either number of segments the cross sections for the
 *     shape are recalculated */
	void setSegments(int x, int y) {
		if (m_segmentsX != x || m_segmentsY != y) {
			m_segmentsX = x;
			m_segmentsY = y;
			if (crossSections != NULL) {
				delete[] crossSections;
				crossSections = NULL;
			}
			createCrossSections();
		}
	}

/* Given two cross sections, draws a connected strip between them */
	void drawConnectedStrip(VertexList *bottom, VertexList *top) {
		if (top->length() == bottom->length()) {
			gVertex *top_v, *bottom_v;
			glBegin(GL_TRIANGLE_STRIP);
			
			for (int i = 0; i < top->length(); ++i) {
				top_v = top->get_Vertex(i);
				glNormal3f(top_v->n[0], top_v->n[1], top_v->n[2]);				
				glVertex3f(top_v->p[0], top_v->p[1], top_v->p[2]);


				bottom_v = bottom->get_Vertex(i);
				glNormal3f(bottom_v->n[0], bottom_v->n[1], bottom_v->n[2]);
				glVertex3f(bottom_v->p[0], bottom_v->p[1], bottom_v->p[2]);

			}

            /* Connects last VertexList to the first VertexList */
			top_v = top->get_Vertex(0);
			glVertex3f(top_v->p[0], top_v->p[1], top_v->p[2]);

			bottom_v = bottom->get_Vertex(0);
			glVertex3f(bottom_v->p[0], bottom_v->p[1], bottom_v->p[2]);

			glEnd();

		} else {
            return;
		}
	}

    /* Given two cross sections, draws a strip between them without connecting
 *     the first and the last VertexList of both. */
	void drawStrip(VertexList *bottom, VertexList *top) {
		if (top->length() == bottom->length()) {
			gVertex *top_v, *bottom_v;
			glBegin(GL_TRIANGLE_STRIP);
			for (int i = 0; i < top->length(); ++i) {

				top_v = top->get_Vertex(i);
				bottom_v = bottom->get_Vertex(i);
				glNormal3f(top_v->n[0], top_v->n[1], top_v->n[2]);
				glVertex3f(top_v->p[0], top_v->p[1], top_v->p[2]);
				glNormal3f(bottom_v->n[0], bottom_v->n[1], bottom_v->n[2]);
				glVertex3f(bottom_v->p[0], bottom_v->p[1], bottom_v->p[2]);

			}
			glEnd();

		} else {
            return;
		}
	}

    /* Given a cross section and a radial point, creates a connected fan
 *     around the point that is connected to each outer vertex */
	void drawConnectedFan(gVertex *center, VertexList *outer) {
		gVertex *v;
		glBegin(GL_TRIANGLE_FAN);
		glNormal3f(center->n[0], center->n[1], center->n[2]);
		glVertex3f(center->p[0], center->p[1], center->p[2]);
		for (int i = 0; i < outer->length(); ++i) {
			v = outer->get_Vertex(i);
			glNormal3f(v->n[0], v->n[1], v->n[2]);
			glVertex3f(v->p[0], v->p[1], v->p[2]);

		}
		v = outer->get_Vertex(0);
		glVertex3f(v->p[0], v->p[1], v->p[2]);
		glEnd();
	}

	virtual void draw() {};
	virtual void drawNormal() {};
	virtual void createCrossSections() {};
	virtual double intersect(Point p, Vector ray) {};
	virtual Vector getIsectNormal(Point p, Vector ray, double t) {};
	virtual Point mapToSquare(Point p) {};

protected:
	void normalizeNormal (float x, float y, float z) {
		normalizeNormal (Vector(x, y, z));
	};

	void normalizeNormal (Vector v) {
		v.normalize();
		glNormal3dv(v.unpack());
	};
    	double get_t(double proj_val, Point p, Vector ray, IntersectPlane plane) {
        	double t = -1.0;
        	switch(plane) {
            	case X:
                	t = (proj_val - p[0]) / (ray[0]);
                	break;
            	case Y:
			t = (proj_val - p[1]) / (ray[1]);
			break;
		    case Z:
			t = (proj_val - p[2]) / (ray[2]);
			break;
		}
		return t;
	    };
	void solve_t(double t, double *x, double *y, double *z, Point p, Vector ray, IntersectPlane plane) {
        	switch(plane) {
        	    case X:
			(*x) = (*x);
			(*y) = p[1] + (ray[1] * t);
			(*z) = p[2] + (ray[2] * t);
			break;
		    case Y:
			(*x) = p[0] + (ray[0] * t);
			(*y) = (*y);
			(*z) = p[2] + (ray[2] * t);
			break;
		    case Z:
			(*x) = p[0] + (ray[0] * t);
			(*y) = p[1] + (ray[1] * t);
			(*z) = (*z);
			break;
		}
    	};



	int m_segmentsX, m_segmentsY;
	double x_dist, y_dist, z_dist;

	int numCrossSections;
	VertexList *crossSections;
};

#endif
