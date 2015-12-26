/* Casey Gowrie & Jake Mingolla
 * Sphere class for assignment 1 
 * Comp 175 -- Remco Chang */

#ifndef SPHERE_H
#define SPHERE_H

#include "Shape.h"

class Sphere : public Shape {
public:
	Sphere() {
		Vector normalTop(0.0, 1.0, 0.0);
		Vector normalBottom(0.0, -1.0, 0.0);
		topCenter = gVertex(0.0, y_dist * .5, 0.0, normalTop);
		bottomCenter = gVertex(0.0, y_dist * -.5, 0.0, normalBottom);
	};
	~Sphere() {};
	
	Vector getIsectNormal(Point p, Vector ray, double t) {
		Point on_surface = p + t * ray;
		//Point origin = Point(0,0,0);
		double x = on_surface[0];
		double y = on_surface[1];
		double z = on_surface[2];
		Vector normal = Vector(x, y, z);
		//normal.normalize();
		return normal;
	};

	Point mapToSquare (Point p)
	{
		Point s = Point(0, 0, 0);
        double x = p[0];
        double y = p[1];
        double z = p[2];

    	double theta = atan2(z, x);

    	s[0] = -theta / (2 * PI) + 0.5;
    	s[1] = -1 * asin(2 * y) / PI + 0.5;

    	return s;
	};
		

	double intersect(Point p, Vector ray) {
		double a, b, c, minus_t, plus_t, t;
		a = ray[0] * ray[0] + ray[1] * ray[1] + ray[2] * ray[2];
		b = 2 * (p[0] * ray[0] + p[1] * ray[1] + p[2] * ray[2]);
		c = p[0] * p[0] + p[1] * p[1] + p[2] * p[2] - .5*.5;
		return solve_quadratic(a, b, c);
	};

	/* draws bottom fan, then strips for each cross sec, then top fan */
	void draw() {
		drawConnectedFan(&bottomCenter, &crossSections[0]);
		for (int i = 0; i < numCrossSections - 1; i++) {
			drawConnectedStrip(&crossSections[i], &crossSections[i+1]);
		}
		drawConnectedFan(&topCenter, &crossSections[numCrossSections - 1]);
	};

    /* Draws the normals for the sphere, starting with the cross sections
 *     and then drawing the top and bottom fans and ceter points. */
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

				/* two points for normal line segment */
				vertex = vertices->get_Vertex(i);
				glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
				glVertex3f(vertex->p[0] + vertex->n[0] * factor,
                           vertex->p[1] + vertex->n[1] * factor, 
                           vertex->p[2] + vertex->n[2] * factor);
				
				glEnd();
			}
		}

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
			
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

    /* Creates the cross sections of the sphere from bottom to top */
	void createCrossSections() {
		double x, y, z, y_offset, theta;
        double y_theta_offset, y_theta;
		double x_radius, z_radius,y_radius, y_r_squared;
		
		numCrossSections = m_segmentsY - 1;
		y_offset = y_dist / m_segmentsY;

        y_theta_offset = PI  / m_segmentsY;
        y_theta = (-0.5 * PI);
        y_radius = (.5 * y_dist);

		Vector normal;

		y_r_squared = (y_radius * y_radius);

		crossSections = new VertexList [numCrossSections];
		for (int c = 0; c < numCrossSections; ++c) {
			y_theta += y_theta_offset;
            y = y_radius * sin(y_theta);
			x_radius = sqrt(y_r_squared - pow(y, 2));
			z_radius = x_radius;
			for (int i = 0; i < m_segmentsX; ++i) {
				theta = ((2 * PI) / m_segmentsX) * i;
				x = x_radius * cos(theta);
				z = z_radius * sin(theta);
				normal = Vector(x, y, z);
				normal.normalize();
				gVertex v(x, y, z, normal);
				
				crossSections[c].add_Vertex(v);
			}
		}
	};
private:
		gVertex topCenter, bottomCenter;
};

#endif
