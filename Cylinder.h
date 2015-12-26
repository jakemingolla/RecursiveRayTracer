/* Jake Mingolla & Casey Gowrie
 * Cylinder class for assignment 1
 * Comp 175 -- Remco Chang */

#ifndef CYLINDER_H
#define CYLINDER_H

#include "Shape.h"
#include <cmath>

class Cylinder : public Shape {
public:
	Cylinder() {
		topFan = NULL;
		bottomFan = NULL;
		Vector normalTop(0.0, 1.0, 0.0);
		Vector normalBottom(0.0, -1.0, 0.0);
		topCenter = gVertex(0.0, y_dist * .5, 0.0, normalTop);
		bottomCenter = gVertex(0.0, y_dist * -.5, 0.0, normalBottom);
	};
	~Cylinder() {
		if (topFan != NULL) {
			delete topFan;
		}
		if (bottomFan != NULL) {
			delete bottomFan;
		}
	};

	Point mapToSquare (Point p)
	{
        Point s = Point(0, 0, 0);
        double x = p[0];
        double y = p[1];
        double z = p[2];

        if (IN_RANGE(y, (double) -0.5) || y == -0.5) {
        	s[0] = -x + 0.5;
            s[1] = z + 0.5;
        } else if (IN_RANGE(y, (double) 0.5) || y == 0.5) {
            s[0] = x + 0.5;
            s[1] = z + 0.5;
        } else {
        	// z = z / .25;
        	// x = x / .25;
        	double theta = atan2(z, x);
        	// double u;
        	// if (theta < 0) {
        	// 	u = -1 * theta / (2 * PI);
        	// } else {
        	// 	u = 1 - (theta / (2 * PI));
        	// }
        	s[0] = -theta / (2 * PI) + 0.5;
        	s[1] = -y + 0.5;
        }
		return s;
	};

	Vector getIsectNormal(Point p, Vector ray, double t) {
		Point on_surface = p + t * ray;
		double x = on_surface[0];
		double y = on_surface[1];
		double z = on_surface[2];

		Vector n = Vector(0, 0, 0);

		if (IN_RANGE_EQ(y, -0.5)) {
			n[1] = -1.0;
		} else if (IN_RANGE_EQ(y, 0.5)) {
			n[1] = 1.0;
		} else {
			n[0] = x;
			n[2] = z;
		}
        n.normalize();
		return n;
	}

	double intersect(Point p, Vector ray) {
		double t, min_t = -1.0;
		double x, y, z;

		y = 0.5;
		t = get_t(y, p, ray, Y);
		solve_t(t, &x, &y, &z, p, ray, Y);
		if (in_cap(x, z)) {
			if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
		}	

		y = -0.5;
		t = get_t(y, p, ray, Y);
		solve_t(t, &x, &y, &z, p, ray, Y);
		if (in_cap(x, z)) {
			if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
		}	

		double A, B, C;
		A = (SQR(ray[2]) + SQR(ray[0]));
		B = 2 * (ray[2] * p[2] + ray[0] * p[0]);
		C = SQR(p[2]) + SQR(p[0]) - SQR(0.5);
		t = solve_quadratic(A, B, C);
		if (!IN_RANGE(t, -1.0)) {
			x = p[0] + (t * ray[0]);
			y = p[1] + (t * ray[1]);
			z = p[2] + (t * ray[2]);
			//if (in_cap(x, z) && in_range(y)) {
			if (in_range(y)) {
				if (t > 0 && (t < min_t || min_t == -1.0))
            		min_t = t;
			}	
		}
		return min_t;
	};

    /* Draws the surface of the cylinder starting from the bottom and
 *     going to the top. */
	void draw() {
		drawConnectedFan(&bottomCenter, bottomFan);
		for (int i = 0; i < numCrossSections - 1; i++) {
			drawConnectedStrip(&crossSections[i], &crossSections[i+1]);
		}
		drawConnectedFan(&topCenter, topFan);
	};

    /* Draws the normals of the cylinder starting at the cross sections
 *     then drawing the top and bottom. */
	void drawNormal() {
		VertexList *vertices;
		gVertex *vertex;
		int length;
		double factor = 0.1f;
		for (int c = 0; c < numCrossSections; ++c) {
			vertices = &crossSections[c];
			length = vertices->length();
			for (int i = 0; i < length; ++i) {
				glBegin(GL_LINES);
				glColor3f(1.0f, 0.0f, 0.0f);
				vertex = vertices->get_Vertex(i);
				glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
				glVertex3f(vertex->p[0] + vertex->n[0] * factor, 
                           vertex->p[1] + vertex->n[1] * factor, 
                           vertex->p[2] + vertex->n[2] * factor);

				glEnd();
			}
		}
		for (int i = 0; i < m_segmentsX; ++i) {
			glBegin(GL_LINES);
			glColor3f(1.0f, 0.0f, 0.0f);
			
			vertex = bottomFan->get_Vertex(i);
			glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
			glVertex3f(vertex->p[0] + vertex->n[0] * factor, 
                       vertex->p[1] + vertex->n[1] * factor, 
                       vertex->p[2] + vertex->n[2] * factor);

			vertex = topFan->get_Vertex(i);
			glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
			glVertex3f(vertex->p[0] + vertex->n[0] * factor,
                       vertex->p[1] + vertex->n[1] * factor, 
                       vertex->p[2] + vertex->n[2] * factor);
			glEnd();
		}

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
			
        /* Also draws the normal for the center points on the top
 *         and the bottom. */
		vertex = &bottomCenter;
		glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
		glVertex3f(vertex->p[0] + vertex->n[0] * factor, 
                   vertex->p[1] + vertex->n[1] * factor, 
                   vertex->p[2] + vertex->n[2] * factor);

		vertex = &topCenter;
		glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
		glVertex3f(vertex->p[0] + vertex->n[0] * factor,
                   vertex->p[1] + vertex->n[1] * factor,
                   vertex->p[2] + vertex->n[2] * factor);
		glEnd();

	};

    /* Creates cross sections and the bottom and top fans of the cylinder */
	void createCrossSections() {
		double x, y, z, y_offset, theta;
        double x_radius = .5 * x_dist, z_radius = .5 * z_dist;
		double max_y = 0.5 * y_dist, min_y = -0.5 * y_dist;
		Vector normal;
		gVertex v;
		numCrossSections = m_segmentsY + 1;
		y_offset = y_dist / m_segmentsY;
		y = min_y;
		crossSections = new VertexList [numCrossSections];
		for (int c = 0; c < numCrossSections; ++c) {
			for (int i = 0; i < m_segmentsX; ++i) {
				theta = ((2 * PI) / m_segmentsX) * i;
				x = x_radius * cos(theta);
				z = z_radius * sin(theta);
				normal = Vector(x, 0.0, z);
				normal.normalize();
				v = gVertex(x, y, z, normal);
				
				crossSections[c].add_Vertex(v);
			}
			y += y_offset;
		}

        if (bottomFan != NULL) {
            delete bottomFan;
        }
		bottomFan = new VertexList;
		for (int i = 0; i < m_segmentsX; ++i) {
			theta = ((2 * PI) / m_segmentsX) * i;
			x = x_radius * cos(theta);
			z = z_radius * sin(theta);
			normal = Vector(0.0, -1.0, 0.0);
			v = gVertex(x, min_y, z, normal);
			bottomFan->add_Vertex(v);
		}
		
        if (topFan != NULL) {
            delete topFan;
        }
		topFan = new VertexList;
		for (int i = 0; i < m_segmentsX; ++i) {
			theta = ((2 * PI) / m_segmentsX) * i;
				x = x_radius * cos(theta);
				z = z_radius * sin(theta);
				normal = Vector(0.0, 1.0, 0.0);
				v = gVertex(x, max_y, z, normal);
				topFan->add_Vertex(v);
		}
	};
private:
	VertexList *topFan, *bottomFan;
	gVertex topCenter, bottomCenter;

	bool in_cap(double x, double z) {
		if((SQR(x) + SQR(z)) <= SQR(0.5)) {
			return true;
		} else {
			return false;
		}
	}
	
	bool in_range(double y) {
		return ((y <= 0.5) && (y >= -0.5));
	}
};
#endif
