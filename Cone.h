/* Jake Mingolla & Casey Gowrie.
 * Cone class for assignment 1
 * Comp 175 -- Remco Chang */

#ifndef CONE_H
#define CONE_H

#include "Shape.h"

class Cone : public Shape {
public:
    /* Default constructor for the Cone class. */
	Cone() {
		bottomFan = NULL;
		Vector normalBottom(0.0, -1.0, 0.0);
		bottomCenter = gVertex(0.0, y_dist * -.5, 0.0, normalBottom);
	};
    /* Destructor for the Cone class. */
	~Cone() {
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
        } else {
        	double theta = atan2(z, x);

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
		} else {
			Vector V = Vector(x, y, z);
			V[1] = 0;
			V.normalize();
			// double m = -2.0;
			// double normal_m = -1 / m;
			// double y_n = 0.5 * normal_m;
			// double scale = SQRT(SQR(x) + SQR(z));

			n[0] = V[0];
			n[1] = .5;
			n[2] = V[2];
			n.normalize();
		}

		return n;
	}

	double intersect(Point p, Vector ray) {
		double t, min_t = -1.0;
		double x, y, z;

		y = -0.5;
		t = get_t(y, p, ray, Y);
		solve_t(t, &x, &y, &z, p, ray, Y);
		if (in_cap(x, z)) {
			if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
		}

		double A, B, C;
		double apex = 0.5;
		// try 1
		// A = SQR(ray[0]) + SQR(ray[2]) + (SQR(0.5) * SQR(ray[1]));
		// B = 2 * ((p[0] * ray[0]) + (p[2] * ray[2]) + (SQR(0.5) * apex * ray[1]) - (SQR(0.5) * ray[1]));
		// C = SQR(p[0]) + SQR(p[2]) - (SQR(0.5) * apex) + (SQR(0.5) * SQR(p[1])) + (2 * SQR(0.5) * apex * p[1]);
	
		// try 2
		// A = SQR(ray[0]) + SQR(ray[2]) - (SQR(0.5) * SQR(ray[1]));
		// B = 2 * ((p[0] * ray[0]) + (p[2] * ray[2]) + SQR(0.5) * ((apex * ray[1]) - (p[1] * ray[1])));
		// C = SQR(p[0]) + SQR(p[1]) + SQR(0.5) * ((-1 * SQR(apex)) + (2 * apex * p[1]) - SQR(p[1]));

		// try 3
		// A = SQR(ray[0]) + SQR(ray[1]) - SQR(ray[2]);
		// B = 2 * (p[0] * ray[0] + p[1] * ray[1] - p[2] * ray[2]);
		// C = SQR(p[0]) + SQR(p[1]) - SQR(p[2]);
		
		// try 4
		A = SQR(ray[0]) + SQR(ray[2]) - (SQR(0.5) * SQR(ray[1]));	
		B = 2 * (p[0] * ray[0] + p[2] * ray[2]) - (SQR(0.5) * 2 * (p[1] * ray[1] - ray[1] * apex));
		C = (SQR(p[0]) + SQR(p[2])) - (SQR(0.5) * (SQR(apex) - 2 * p[1] * apex + SQR(p[1])));
		

		t = solve_quadratic(A, B, C);
		if (!IN_RANGE(t, -1.0)) {
			x = p[0] + (t * ray[0]);
			y = p[1] + (t * ray[1]);
			z = p[2] + (t * ray[2]);
			if (in_range(y)) {
				if (t > 0 && (t < min_t || min_t == -1.0))
            		min_t = t;
			}	
		}
	

		return min_t;	
	};


    /* Draws the cone by drawing connected strip between each crossSection then
 *     drawing a connected fan for the bottom base */
	void draw() {
		drawConnectedFan(&bottomCenter, bottomFan);
		for (int i = 0; i < numCrossSections - 1; i++) {
			drawConnectedStrip(&crossSections[i], &crossSections[i+1]);
		}
	};

    /* Draws the normal for the Cone, iterating first through the crossSection
 *     VertexLists then by drawing the normals of the bottom. */
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

			glEnd();
		}

		glBegin(GL_LINES);
		glColor3f(1.0f, 0.0f, 0.0f);
			
		vertex = &bottomCenter;
		glVertex3f(vertex->p[0], vertex->p[1], vertex->p[2]);
		glVertex3f(vertex->p[0] + vertex->n[0] * factor, 
                   vertex->p[1] + vertex->n[1] * factor, 
                   vertex->p[2] + vertex->n[2] * factor);
		glEnd();
	};

    /* Creates the cross sections for the Cone, populating both the
 *     crossSections array of VertexList as well as the bottom fan
 *     VertexList. */
	void createCrossSections() {
		double x, y, z, y_offset, theta, x_base_radius = .5 * x_dist;
        /* m is the slope of the line creating the outer edge of the cone */
        double min_y = -0.5 * y_dist, m = -2.0 * (y_dist / x_dist);
        double z_base_radius = .5 * z_dist, normal_m = -1.0 / m;
		double height_from_base, height_from_top, x_radius, z_radius;
   
		Vector normal;
		gVertex v;
        double y_n = normal_m * x_base_radius;
		numCrossSections = m_segmentsY + 1;
		y_offset = y_dist / m_segmentsY;
		
        /* Start creating the cross sections from the base. */
		y = min_y;
		crossSections = new VertexList [numCrossSections];
		for (int c = 0; c < numCrossSections; ++c) {
			height_from_base = c * y_offset;
			height_from_top = y_dist - height_from_base;
			x_radius = (x_base_radius / y_dist) * height_from_top;
			z_radius = (z_base_radius / y_dist) * height_from_top;
			for (int i = 0; i < m_segmentsX; ++i) {
				theta = ((2 * PI) / m_segmentsX) * i;
				x = x_radius * cos(theta);
				z = z_radius * sin(theta);


                normal = Vector(x_base_radius * cos(theta), 
                                y_n,
                                z_base_radius * sin(theta));
                normal.normalize();
    
				v = gVertex(x, y, z, normal);
				crossSections[c].add_Vertex(v);
			}
			y += y_offset;
		}

        /* Recreate bottom fan for new m_segmentsX */
		if (bottomFan != NULL) {
			delete bottomFan;
		}
		bottomFan = new VertexList;
		for (int i = 0; i < m_segmentsX; ++i) {
			theta = ((2 * PI) / m_segmentsX) * i; 
			x = x_base_radius * cos(theta);
			z = z_base_radius * sin(theta);
            /* all normals on the bottom point directly down. */
			normal = Vector(0.0, -1.0, 0.0);
			v = gVertex(x, min_y, z, normal);
			bottomFan->add_Vertex(v);
		}
	};
private:
	VertexList *bottomFan;
	gVertex bottomCenter;
	
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
