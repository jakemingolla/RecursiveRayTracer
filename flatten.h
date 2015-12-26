#ifndef FLATTEN_H
#define FLATTEN_H

#include "Algebra.h"
#include "SceneData.h"
#include <vector>
#include "ppm.h"

struct FlatInfo {
	Matrix transformation;
	std::vector<ScenePrimitive*> primitives;
	int intersectIndex;
};

struct FlatNode {
	FlatInfo info;
	FlatNode *next;
};

struct TexInfo {
	string filename;
    int width, height, len;
    char *arr;
};

struct TexNode {
	TexInfo info;
	TexNode *next;
};

FlatNode *flatRoot = NULL;
TexNode *texRoot = NULL;
string currFilenamePath = " ";

void printMatrix1(Matrix m) {
	m = transpose(m);
	std::cout << m[0] << " " << m[1] << " " << m[2] << " " << m[3] << std::endl;
	std::cout << m[4] << " " << m[5] << " " << m[6] << " " << m[7] << std::endl;
	std::cout << m[8] << " " << m[9] << " " << m[10] << " " << m[11] << std::endl;
	std::cout << m[12] << " " << m[13] << " " << m[14] << " " << m[15] << std::endl;
	std::cout << std::endl;
}

void printList(FlatNode *curr) {
	if (curr == NULL)
		return;
	printMatrix1((curr->info).transformation);
	printList(curr->next);

}


Matrix makeTranformM(SceneTransformation *t) {
	Matrix M;
	switch(t->type) {
		case TRANSFORMATION_TRANSLATE:
			M = trans_mat(t->translate);
			break;
		case TRANSFORMATION_ROTATE:
			M = rot_mat(t->rotate, t->angle);
			break;
		case TRANSFORMATION_SCALE:
			M = scale_mat(t->scale);
			break;
		case TRANSFORMATION_MATRIX:
			M = t->matrix;
			break;
	}

	return M;
}

void flatten(SceneNode *head, Matrix parentTransform) {
	//Matrix updated;
	//base case doesn't add anything to the global linked list, might need some change
	if (head == NULL)
		return;

	FlatInfo currInfo;
	FlatNode *currNode = new FlatNode;

	//store all primitives, as attempting to do below, check on vector info
	int len = head->primitives.size();
	currInfo.primitives = std::vector<ScenePrimitive*>(len);
	for (int i = 0; i < len; ++i) {
		currInfo.primitives[i] = new ScenePrimitive;
		(currInfo.primitives)[i]->type = (head->primitives[i])->type;
		(currInfo.primitives)[i]->meshfile = (head->primitives[i])->meshfile;
		(currInfo.primitives)[i]->material = (head->primitives[i])->material;
	}

	//iterate through transformations at that node
	//if a vector, use info in vector to make curr's matrix (Algebra.h, see camera rotate things for example)
	len = head->transformations.size();
	for (int i = 0; i < len; ++i) {
		currInfo.transformation = currInfo.transformation * makeTranformM(head->transformations[i]);
	}

	currInfo.transformation = parentTransform * currInfo.transformation;

	//curr Tranform matrix = parent * one created above (product of transformations looped through)
	currNode->info = currInfo;
	currNode->next = flatRoot;
	flatRoot = currNode;
	

	//loop through children of SceneNode, recurse on all, with curr Tranform matrix as parentTransform arg
	len = head->children.size();
	for (int i = 0; i < len; ++i) {
		flatten(head->children[i], currInfo.transformation);
	}
}

void deleteList() {
	FlatNode *iter = flatRoot;
	FlatNode *prev = NULL;
	while (iter != NULL) {
		prev = iter;
		iter = iter->next;
		int len = (prev->info).primitives.size();
		for (int i = 0; i < len; ++i) {
			delete prev->info.primitives[i];
		}
		delete prev;
	}
	flatRoot = NULL;
}

void flatten(SceneNode *head, bool reflatten) {
	if (reflatten) {
		if (head != NULL) {
			deleteList();
			flatRoot = NULL;
		}
		flatten(head, Matrix());
	}
}

void renderMesh(string file) {
	(void) file;
	return;
}

TexInfo *findTexture(string filename) {
	TexNode *iter = texRoot;
	while (iter != NULL) {
		TexInfo *info = &(iter->info);
		if (info->filename == filename) {
//                 std::cout << "filename = " << info->filename << endl;
//                 char *tex_array = info->arr;
//                 ppm test = ppm("./data/tests/image/earth.ppm");
// 
//                 std::cout << "dim = " << test.getWidth() * test.getHeight() * 3  << " == " << info->len << endl;
//                 char *arr = test.getPixels();
//         
//                 for (int i = 0; i < test.getHeight() * test.getWidth() * 3; ++i) {
//                     if (arr[i] != tex_array[i]) {
//                         std::cout << "not equal! arr[" << i << "] = " << (int)arr[i] << ", tex_array[" << i << "] = " << (int)tex_array[i] << endl;
//                     }
//                     else { cout << i << " okay" << endl; }
//                 }
//                 std::cout << "done!" << endl;
//                 exit(1);




			return info;
		}
		iter = iter->next;
	}
	return NULL;
}

void deleteTextureList() {
	TexNode *iter = texRoot;
	TexNode *prev = NULL;
	while (iter != NULL) {
		prev = iter;
		iter = iter->next;
		delete[] prev->info.arr;
		delete prev;
	}
	texRoot = NULL;
}

void flattenTextures() {
	FlatNode *iter = flatRoot;
	while (iter != NULL) {
		FlatInfo info = iter->info;
		int len = info.primitives.size();
		for (int i = 0; i < len; ++i) {
			SceneMaterial mat = info.primitives[i]->material;
			SceneFileMap *tex = mat.textureMap;
			if ((tex->isUsed) && (findTexture(tex->filename) == NULL)) {
				TexNode *new_tex = new TexNode;

			 	ppm tex_ppm = ppm(tex->filename);
                char *colors = tex_ppm.getPixels();
                int len = tex_ppm.getWidth() * tex_ppm.getHeight() * 3;
                new_tex->info.arr = new char[len];

                for (int i = 0; i < len; ++i) {
                    new_tex->info.arr[i] = colors[i];
                }

                new_tex->info.len = len;
                new_tex->info.height = tex_ppm.getHeight();
                new_tex->info.width = tex_ppm.getWidth();
			 	new_tex->info.filename = tex->filename;
			 	new_tex->next = texRoot;
			 	texRoot = new_tex;
			 }
		}
		iter = iter->next;
	}
}



#endif
