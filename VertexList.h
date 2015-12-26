/* Casey Gowrie, Jake Mingolla
 * gVertex.h
 */

#ifndef VERTEXLIST_H
#define VERTEXLIST_H

 #include "Algebra.h"

/* struct for a vertex, whcih contains a point and a normal vector */
 struct gVertex {
	Point p;
 	Vector n;

 	gVertex() {
 		p = Point();
 		n = Vector();
 	}

 	gVertex(double init_x, double init_y, double init_z) {
 		p = Point(init_x, init_y, init_z);
 		n = Vector(0.0f, 0.0f, 0.0f);
 	}

 	gVertex(double init_x, double init_y, double init_z, Vector v) {
 		p = Point(init_x, init_y, init_z);
 		n = v;
 	}
 };


/*vertex list class
 allows for basic operations on a list of vertices */
 class VertexList {
 	public:
 		VertexList() {
 			count = 0;
 			capacity = 8;
 			arr = new gVertex [capacity];
 		}

 		~VertexList() {
 			if (arr != NULL) {
  				delete[] arr;
  			}
 		}

 		int add_Vertex(gVertex v) {
 			int willSucceed = 1;
 			if (count >= capacity) {
 				willSucceed = expand();
 			}
 			if (willSucceed) {
 				arr[count] = v;
 				++count;
 			}
 			return willSucceed;
 		}

 		int length() { return count; }

 		gVertex *get_Vertex(int i) {
 			return (i >= 0 && i < count) ? &arr[i] : NULL;
 		}



 	private:
 		int capacity;
 		int count;
 		gVertex *arr;

 		int expand() {
 			gVertex *temp = NULL;

 			temp = new gVertex [capacity * 2];

 			for (int i = 0; i < count; ++i) {
 				temp[i] = arr[i];
 			}

 			delete[] arr;
 			arr = temp;
 			capacity *= 2;

 			return 1;
 		}
 };

 #endif
