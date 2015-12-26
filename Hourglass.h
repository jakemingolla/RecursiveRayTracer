/* Jake Mingolla & Casey Gowrie
 * Special Shape 1 (Hourglass) class for assignment 1
 * Comp 175 -- Remco Chang */
#ifndef HOURGLASS_H
#define HOURGLASS_H

#include "Shape.h"

class Hourglass : public Shape {
public:
	Hourglass() {
		topFan = NULL;
		bottomFan = NULL;
		Vector normalTop(0.0, 1.0, 0.0);
		Vector normalBottom(0.0, -1.0, 0.0);
		topCenter = gVertex(0.0, y_dist * .5, 0.0, normalTop);
		bottomCenter = gVertex(0.0, y_dist * -.5, 0.0, normalBottom);
	};
	~Hourglass() {
		if (topFan != NULL) {
			delete topFan;
		}
		if (bottomFan != NULL) {
			delete bottomFan;
		}
        if (crossSections != NULL) {
            delete[] crossSections;
            crossSections = NULL;
        }
    };

    /* Draws the hourglass from bottom to top */
	void draw() {
		drawConnectedFan(&bottomCenter, bottomFan);
		for (int i = 0; i < numCrossSections - 1; i++) {
			drawConnectedStrip(&crossSections[i], &crossSections[i+1]);
		}
		drawConnectedFan(&topCenter, topFan);
	};

    /* Draws the normals, first of the cross sections then of the bottom
 *     and top fans. */
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
			
        /* Also draws the normals for the bottom and top center points */
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

	void createCrossSections() {
		double x, y, z, y_offset, theta;
        double x_radius = .5 * x_dist, z_radius = .5 * z_dist;
		double y_r_squared, sphere_radius;
        double max_y = 0.5 * y_dist, min_y = -0.5 * y_dist;
		Vector normal;
		gVertex v;
		numCrossSections = m_segmentsY + 1;
		y_offset = y_dist / m_segmentsY;
		y = min_y;
        /* y radius squared is used to find the sphere radius at each
 *         step along the hourglass */
		y_r_squared = pow(.5 * y_dist, 2);
		crossSections = new VertexList [numCrossSections];
		for (int c = 0; c < numCrossSections - 1; ++c) {
			sphere_radius = sqrt(y_r_squared - (y * y));
			x_radius = .5 * x_dist - sphere_radius;
			z_radius = .5 * z_dist - sphere_radius;
			for (int i = 0; i < m_segmentsX; ++i) {
				theta = ((2 * PI) / m_segmentsX) * i;
				x = x_radius * cos(theta);
				z = z_radius * sin(theta);
				normal = Vector((.5 * x_dist) * cos(theta),
                                 -1 * y,
                                (.5 * z_dist) * sin(theta));
				normal.normalize();
				v = gVertex(x, y, z, normal);
				
				crossSections[c].add_Vertex(v);
			}
			y += y_offset;
		}

    /* Creates the last cross section to avoid a bug in which the top cross
 *     section of the hourglass was only displayed with msegments_Y was a power
 *     of 2. While this a loss of efficiency, it works. */
		x_radius = .5 * x_dist;
		z_radius = .5 * z_dist;
		for (int i = 0; i < m_segmentsX; ++i) {
			theta = ((2 * PI) / m_segmentsX) * i;
			x = x_radius * cos(theta);
			z = z_radius * sin(theta);
			normal = Vector((.5 * x_dist) * cos(theta), 
                            -1 * y,
                            (.5 * z_dist) * sin(theta));
			normal.normalize();
			v = gVertex(x, y, z, normal);
			
			crossSections[numCrossSections - 1].add_Vertex(v);
		}

		x_radius = .5 * x_dist;
		z_radius = .5 * z_dist;

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
};
#endif
