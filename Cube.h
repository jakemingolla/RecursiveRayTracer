/* Jake Mingolla & Casey Gowrie
 * Cube class for assignment 1
 * Comp 175 -- Remco Chang */
#ifndef CUBE_H
#define CUBE_H

#include "Shape.h"

class Cube : public Shape {
public:
    /* Default constructor for the cube class */
	Cube() {
		topStrips = NULL;
		bottomStrips = NULL;
        numStrips = 0;
	};
    /* Destructor for the cube class */
	~Cube() {
		if (topStrips != NULL) {
			delete[] topStrips;
		}
		if (bottomStrips != NULL) {
			delete[] bottomStrips;
		}
	};

	Point mapToSquare(Point p) {
		Point s = Point(0, 0, 0);

		double x = p[0];
		double y = p[1];
		double z = p[2];


		if (IN_RANGE(x, (double) -0.5) || x == -0.5) {
			s[0] = -z;
			s[1] = -y;
		} else if (IN_RANGE(x, (double) 0.5) || x == 0.5) {
			s[0] = -z;
			s[1] = -y;
		} else if (IN_RANGE(y, (double) -0.5) || y == -0.5) {
			s[0] = x;
			s[1] = -z;
		} else if (IN_RANGE(y, (double) 0.5) || y == 0.5) {
			s[0] = x;
			s[1] = -z;
		} else if (IN_RANGE(z, (double) -0.5) || z == -0.5) {
			s[0] = x;
			s[1] = -y;
		} else if (IN_RANGE(z, (double) 0.5) || z == 0.5) {
			s[0] = x;
			s[1] = -y;
		} else {
			std::cerr << "Not on cube: (" << x << ", " << y  << ", " << z << ")" << std::endl;
		}
		// convert to 0-1
		s[0] += 0.5;
		s[1] += 0.5;

		return s;
	};

	Vector getIsectNormal(Point p, Vector ray, double t) {
		Point on_surface = p + t * ray;
		double x = on_surface[0];
		double y = on_surface[1];
		double z = on_surface[2];
		Vector n = Vector(0, 0, 0);
		if (IN_RANGE(x, (double) -0.5) || x == -0.5) {
			n[0] = -1.0;
		} else if (IN_RANGE(x, (double) 0.5) || x == 0.5) {
			n[0] = 1.0;
		} else if (IN_RANGE(y, (double) -0.5) || y == -0.5) {
			n[1] = -1.0;
		} else if (IN_RANGE(y, (double) 0.5) || y == 0.5) {
			n[1] = 1.0;
		} else if (IN_RANGE(z, (double) -0.5) || z == -0.5) {
			n[2] = -1.0;
		} else if (IN_RANGE(z, (double) 0.5) || z == 0.5) {
			n[2] = 1.0;
		} else {
			std::cerr << "No normal found for: (" << x << ", " << y  << ", " << z << ")" << std::endl;
		}
		return n;
	};

    double intersect(Point p, Vector ray) {
        double t, x, y, z;
        double min_t = -1.0;

        x = -0.5;
        t = get_t(x, p, ray, X);
        solve_t(t, &x, &y, &z, p, ray, X);
        if (in_bounds(y) && in_bounds(z)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 


        x = 0.5;
        t = get_t(x, p, ray, X);
        solve_t(t, &x, &y, &z, p, ray, X);
        if (in_bounds(y) && in_bounds(z)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 


        y = -0.5;
        t = get_t(y, p, ray, Y);
        solve_t(t, &x, &y, &z, p, ray, Y);
        if (in_bounds(x) && in_bounds(z)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 

        y = 0.5;
        t = get_t(y, p, ray, Y);
        solve_t(t, &x, &y, &z, p, ray, Y);
        if (in_bounds(x) && in_bounds(z)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 


        z = -0.5;
        t = get_t(z, p, ray, Z);
        solve_t(t, &x, &y, &z, p, ray, Z);
        if (in_bounds(x) && in_bounds(y)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 


        z = 0.5;
        t = get_t(z, p, ray, Z);
        solve_t(t, &x, &y, &z, p, ray, Z);
        if (in_bounds(x) && in_bounds(y)) {
            if (t > 0 && (t < min_t || min_t == -1.0))
            	min_t = t;
        } 

	// no intersect
        return min_t;
    };

    

    /* Draws the cube by first drawing disconnected strips for the top and
 *     the bottom of the cube then by drawing connected strips for each
 *     cross section */
	void draw() {
		for (int i = 0; i < numStrips - 1; ++i) {
			drawStrip(&topStrips[i], &topStrips[i+1]);
			drawStrip(&bottomStrips[i], &bottomStrips[i+1]);
		}
		for (int i = 0; i < numCrossSections - 1; ++i) {
			drawConnectedStrip(&crossSections[i], &crossSections[i+1]);
		}
	};


    /* Draws the cube normals by first drawing disconnected strips for the top
 *     and the bottom of the cube then by drawing connected strips for each
 *     cross section */
	void drawNormal() {
		VertexList *vertices;
		gVertex *vertex;
		int length;
		double factor = 0.1f;
		for (int c = 0; c < numStrips; ++c) {
			vertices = &topStrips[c];
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
			vertices = &bottomStrips[c];
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

	};

/* Creates the cross sections of the cube by designating points along
 * the x-z plane. Looking down on the x-z plane along the square
 * (-.5, -.5) to (.5, .5), it creates points by going around
 * from the bottom left. */
	void createCrossSections() {
		double x, y, z, x_offset, y_offset, z_offset;
        double max_y = y_dist / 2, min_y = -1 * max_y;
		double max_x = .5 * x_dist, min_x = -1 * max_x;
        double  max_z = .5 * z_dist, min_z = -1 * max_z;

		Vector normal;

		numStrips = m_segmentsX + 1;

		if (topStrips != NULL) {
			delete[] topStrips;
		}
		topStrips = new VertexList [numStrips];
		
		if (bottomStrips != NULL) {
			delete[] bottomStrips;
		}
		bottomStrips = new VertexList [numStrips];

		x_offset = x_dist / m_segmentsX;
		y_offset = y_dist / m_segmentsY;
		z_offset = z_dist / m_segmentsX;


		x = min_x;
		z = min_z;
        /* Populates the bottom and the top */
		for (int i = 0; i < numStrips; ++i) {
			for (int k = 0; k < numStrips; ++k) {
				gVertex top(x, max_y, z, Vector(0.0, 1.0, 0.0));
				gVertex bottom(x, min_y, z, Vector(0.0, -1.0, 0.0));
				topStrips[i].add_Vertex(top);
				bottomStrips[i].add_Vertex(bottom);
				z += z_offset;
			}
			z = min_z;
			x += x_offset;
		}

		x = min_x;
		y = min_y;
		z = min_z;
		numCrossSections = m_segmentsY + 1;
		crossSections = new VertexList [numCrossSections];
		gVertex v;
		for (int c = 0; c < numCrossSections; ++c) {
            /* Goes along the bottom side of the x-z plane */
			normal = Vector(-1.0, 0.0, 0.0);
			for (int i = 0; i < m_segmentsX + 1; ++i) {
				v.p = Point(min_x, y, z);
				v.n = normal;
				crossSections[c].add_Vertex(v);
				z += z_offset;
			}
			z -= z_offset;

            /* Goes along the right side of the x-z plane */
			normal = Vector(0.0, 0.0, 1.0);
			for (int k = 0; k < m_segmentsX + 1; ++k) {
				v.p = Point(x, y, max_z);
				v.n = normal;
				crossSections[c].add_Vertex(v);
				x += x_offset;
			}
			x -= x_offset;

            /* Goes along the top side of the x-z plane */
			normal = Vector(1.0, 0.0, 0.0);
			for (int i = 0; i < m_segmentsX + 1; ++i) {
				v.p = Point(max_x, y, z);
				v.n = normal;
				crossSections[c].add_Vertex(v);
				z -= z_offset;
			}
			z += z_offset;

            /* Goes along the left side of the x-z plane */
			normal = Vector(0.0, 0.0, -1.0);
			for (int k = 0; k < m_segmentsX + 1; ++k) {
				v.p = Point(x, y, min_z);
				v.n = normal;
				crossSections[c].add_Vertex(v);
				x -= x_offset;
			}
			x += x_offset;
            /* Advances up to the next cross section */
			y += y_offset;
		}


	};
private:
	int numStrips;
	VertexList *topStrips;
	VertexList *bottomStrips;
    
        bool in_bounds(double val) {
        // greater than or equal to?
        if (val <= 0.5 && val >= -0.5) {
            return true;
        } else {
            return false;
        }
    };

};

#endif
