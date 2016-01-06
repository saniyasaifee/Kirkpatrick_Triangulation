/*
 ============================================================================
 Name        : triangulation.c
 Author      : Saniya Saifee
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>

typedef struct point_loc{
	int** triangle;
	int** points;
	struct point_loc* next;
	struct point_loc* previous;
	int tri_count;
}ploc_t;

// A structure to represent an adjacency list node
typedef struct AdjListNode{
	int triangle_no;
    int visited; //used when polygon is created
	struct AdjListNode* next;
}adj_list_node_t;

// A structure to represent an adjacency list
typedef struct AdjList{
	int pt_degree;
	int in_indep_set;
	adj_list_node_t* head;
}adj_list_t;

// A structure to represent a graph. A graph is an array of adjacency lists.
// Size of array will be V (number of vertices in graph)
typedef struct Graph{
	int V;
	adj_list_t* array;
}graph_t;

typedef struct Polygon{
	int num_vertices;
	int* vertex_array;
}polygon_t;
adj_list_node_t* newAdjListNode(int triangle_no);
polygon_t* newPolygon(int poly_vertices);
graph_t* createGraph(int V);
void addPointTriangle(graph_t* graph, int point, int triangle);
int getOrientation(int p1_x, int p1_y, int p2_x, int p2_y, int p3_x, int p3_y);
int find_min_X(polygon_t* polygon,int (*points)[2]);
int is_point_in_triangle(int orient1, int orient2, int orient3 );
int trianngulate_3(polygon_t* polygon,int (*points)[2],int** triangles,graph_t* graph, int ind_point,int new_triangle_count);
int trianngulate_4(polygon_t* polygon, int (*points)[2], int** triangles,graph_t* graph, int ind_point,int new_triangle_count);
int trianngulate_5(polygon_t* polygon, int (*points)[2],int** triangles, graph_t* graph, int ind_point,int new_triangle_count);
int trianngulate_6(polygon_t* polygon, int (*points)[2], int** triangles, graph_t* graph, int ind_point,int new_triangle_count);
ploc_t* create_ploc(int (*points)[2],int (*triangles)[3], int n, int m);
int query_ploc(ploc_t* pl,int x, int y);

// A utility function to create a new adjacency list node
adj_list_node_t* newAdjListNode(int triangle_no){
	adj_list_node_t* newNode = (adj_list_node_t*) malloc (sizeof(adj_list_node_t));
	newNode->triangle_no = triangle_no;
    newNode->visited = 0;
	newNode->next = NULL;
	return newNode;
}
polygon_t* newPolygon(int poly_vertices){
	polygon_t* newPoly = (polygon_t*) malloc (sizeof(polygon_t));
	newPoly->num_vertices = poly_vertices;
	newPoly->vertex_array = (int*)malloc(poly_vertices * sizeof(int));
	return newPoly;
}

// A utility function that creates a graph of V vertices
graph_t* createGraph(int V){
	graph_t* graph = (graph_t*) malloc(sizeof(graph_t));
	graph->V = V;
	// Create an array of adjacency lists.  Size of array will be V
	graph->array = (adj_list_t*)malloc(V * sizeof(adj_list_t));
	// Initialize each adjacency list as empty by making head as NULL
	int i =0;
	for(i=0; i<V; i++){
		graph->array[i].head = NULL;
		graph->array[i].pt_degree = 0;
		graph->array[i].in_indep_set = -1;
	}
	return graph;
}

void addPointTriangle(graph_t* graph, int point, int triangle){
	// Addthe triangle that contains the point.  A new node is added to the adjacency
    // list of point.  The node is added at the begining
	adj_list_node_t* newNode = newAdjListNode(triangle);
	newNode->next = graph->array[point].head;
	graph->array[point].head = newNode;
}
int getOrientation(int p1_x, int p1_y, int p2_x, int p2_y, int p3_x, int p3_y){
	int result_f=0;
	int result = p1_x * (p2_y - p3_y) + p2_x * (p3_y - p1_y) + p3_x * (p1_y - p2_y);
	if (result > 0) {
		result_f = 1;
	} else if (result < 0) {
		result_f = -1;
	} else if (result == 0) {
		result_f = 0;
	}

	return result_f;
}

int find_min_X(polygon_t* polygon,int (*points)[2]){
	int i, point_index,x;
	int min_x_index = 0,current_min = INT_MAX;
	for(i=0; i< (polygon->num_vertices); i++){
		point_index = polygon->vertex_array[i];
		x = points[point_index][0];
		if(x < current_min){
			current_min = x;
			min_x_index = i;
		}
	}
	return min_x_index;
}


int is_point_in_triangle(int orient1, int orient2, int orient3){
	int flag = 0;
	if ((orient1 == orient2) && (orient2==orient3)) {
		flag = 1;
	} else if (orient1 == 0 || orient3 == 0 || orient2 == 0) {

		if (orient2 == 0) {
			if (orient1 == orient3) {
				flag = 1;
			} else {
				flag = 0;
			}
		} else if (orient1 == 0) {
			if (orient2 == orient3) {
				flag = 1;
			} else {
				flag = 0;
			}
		} else if (orient3 == 0) {
			if (orient1 == orient2) {
				flag = 1;
			} else {
				flag = 0;
			}
		}

	} else {
		flag = 0;
	}

	return flag;
}


int trianngulate_3(polygon_t* polygon,int (*points)[2],int** triangles,graph_t* graph, int ind_point,int new_triangle_count){
	int num_vertices,i=0;
	adj_list_node_t *pCrawl;
	num_vertices = polygon->num_vertices;
	int determinant =0;
	int a,b,c;
	if(num_vertices !=3){
		printf("triangulate_3() retriangulate polygon with 3 vertices\n");
		return -1;
	}
	else{
		a=(polygon->vertex_array)[0];
		b=(polygon->vertex_array)[1];
		c=(polygon->vertex_array)[2];
		determinant = getOrientation(points[a][0],points[a][1],points[b][0],points[b][1],points[c][0],points[c][1]);
		triangles[new_triangle_count][i++] = a;
		triangles[new_triangle_count][i++] = b;
		triangles[new_triangle_count][i++] = c;
		pCrawl = graph->array[ind_point].head;
		while(pCrawl){
			triangles[new_triangle_count][i++] = pCrawl->triangle_no;
			pCrawl = pCrawl->next;
		}
	}
	new_triangle_count++;
	return new_triangle_count;
}
int trianngulate_4(polygon_t* polygon, int (*points)[2], int** triangles,graph_t* graph, int ind_point,int new_triangle_count){
	int min_x_index = 0;//pi
	int orient1, orient2, orient3;
	int is_pt_in_tri;
	int i;
	polygon_t* pl_tri_1, *pl_tri_2;
	polygon_t* re_poly;
	re_poly = newPolygon(polygon->num_vertices);
	min_x_index = find_min_X(polygon,points);
	for(i=0; i<polygon->num_vertices; i++){
		re_poly->vertex_array[i] = polygon->vertex_array[min_x_index];
		min_x_index++;
		if(min_x_index == polygon->num_vertices){
			min_x_index = 0;
		}
	}
	int orient_repoly = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1] );
	if(orient_repoly <0){
		int x = re_poly->vertex_array[1];
		int y = re_poly->vertex_array[2];
		int z= re_poly->vertex_array[3];
		re_poly->vertex_array[1] = z;
		re_poly->vertex_array[2] = y;
		re_poly->vertex_array[3] = x;
	}
	pl_tri_1 = newPolygon(3);
	pl_tri_2 = newPolygon(3);
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1] );
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1] );
	orient3 = getOrientation(points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1] );
	is_pt_in_tri = is_point_in_triangle(orient1,orient2, orient3);
	if(is_pt_in_tri == 1){
		//point is in triangle
		pl_tri_1->vertex_array[0] = re_poly->vertex_array[0];
		pl_tri_1->vertex_array[1] = re_poly->vertex_array[2];
		pl_tri_1->vertex_array[2] = re_poly->vertex_array[1];
		pl_tri_2->vertex_array[0] = re_poly->vertex_array[0];
		pl_tri_2->vertex_array[1] = re_poly->vertex_array[2];
		pl_tri_2->vertex_array[2] = re_poly->vertex_array[3];

	}
	else{
		// point is not in triangle
		pl_tri_1->vertex_array[0] = re_poly->vertex_array[0];
		pl_tri_1->vertex_array[1] = re_poly->vertex_array[1];
		pl_tri_1->vertex_array[2] = re_poly->vertex_array[3];
		pl_tri_2->vertex_array[0] = re_poly->vertex_array[1];
		pl_tri_2->vertex_array[1] = re_poly->vertex_array[2];
		pl_tri_2->vertex_array[2] = re_poly->vertex_array[3];
	}

	new_triangle_count = trianngulate_3(pl_tri_1, points,triangles, graph, ind_point, new_triangle_count);
	new_triangle_count = trianngulate_3(pl_tri_2, points,triangles, graph, ind_point, new_triangle_count);
	return new_triangle_count;
}
int trianngulate_5(polygon_t* polygon, int (*points)[2], int** triangles, graph_t* graph, int ind_point,int new_triangle_count){
	polygon_t* re_poly;
	polygon_t* tri_poly;
	polygon_t* fourgon;
	int p1_x,p2_x, p1_y, p2_y;
	int orient1, orient2, orient3;
	int i=0;
	int min_x_index = 0;//pi
	int is_p1_in_tri, is_p2_in_tri;
	min_x_index = find_min_X(polygon,points);
	tri_poly = newPolygon(3);
	fourgon = newPolygon(4);
	re_poly = newPolygon(polygon->num_vertices);
	//reshuffle points in polygon array
	for(i=0; i<polygon->num_vertices; i++){
		re_poly->vertex_array[i] = polygon->vertex_array[min_x_index];
		min_x_index++;
		if(min_x_index == polygon->num_vertices){
			min_x_index = 0;
		}
	}
	int orient_repoly = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1] );
	if(orient_repoly <0){
		int w = re_poly->vertex_array[1];
		int x = re_poly->vertex_array[2];
		int y = re_poly->vertex_array[3];
		int z= re_poly->vertex_array[4];
		re_poly->vertex_array[1] = z;
		re_poly->vertex_array[2] = y;
		re_poly->vertex_array[3] = x;
		re_poly->vertex_array[4] = w;
	}
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	orient3 = getOrientation(points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	is_p1_in_tri = is_point_in_triangle(orient1, orient2, orient3);
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	orient3 = getOrientation(points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	is_p2_in_tri = is_point_in_triangle(orient1, orient2, orient3);
	p1_x = points[re_poly->vertex_array[2]][0];
	p2_x = points[re_poly->vertex_array[3]][0];
	p1_y = points[re_poly->vertex_array[2]][1];
	p2_y = points[re_poly->vertex_array[3]][1];
	if(is_p1_in_tri == 1 && is_p2_in_tri == 1){
		if(p1_x < p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[1];
			tri_poly->vertex_array[2] = re_poly->vertex_array[2];
			fourgon->vertex_array[0] = re_poly->vertex_array[0];
			fourgon->vertex_array[1] = re_poly->vertex_array[2];
			fourgon->vertex_array[2] = re_poly->vertex_array[3];
			fourgon->vertex_array[3] = re_poly->vertex_array[4];
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p1_x > p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[3];
			tri_poly->vertex_array[2] = re_poly->vertex_array[4];
			fourgon->vertex_array[0] = re_poly->vertex_array[0];
			fourgon->vertex_array[1] = re_poly->vertex_array[1];
			fourgon->vertex_array[2] = re_poly->vertex_array[2];
			fourgon->vertex_array[3] = re_poly->vertex_array[3];
			new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p1_x == p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[1];
			tri_poly->vertex_array[2] = re_poly->vertex_array[2];
			fourgon->vertex_array[0] = re_poly->vertex_array[0];
			fourgon->vertex_array[1] = re_poly->vertex_array[2];
			fourgon->vertex_array[2] = re_poly->vertex_array[3];
			fourgon->vertex_array[3] = re_poly->vertex_array[4];
			new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
	}
	else if(is_p1_in_tri == 1 && is_p2_in_tri == 0){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[1];
		tri_poly->vertex_array[2] = re_poly->vertex_array[2];
		fourgon->vertex_array[0] = re_poly->vertex_array[0];
		fourgon->vertex_array[1] = re_poly->vertex_array[2];
		fourgon->vertex_array[2] = re_poly->vertex_array[3];
		fourgon->vertex_array[3] = re_poly->vertex_array[4];
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if(is_p1_in_tri == 0 && is_p2_in_tri == 1){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[3];
		tri_poly->vertex_array[2] = re_poly->vertex_array[4];
		fourgon->vertex_array[0] = re_poly->vertex_array[0];
		fourgon->vertex_array[1] = re_poly->vertex_array[1];
		fourgon->vertex_array[2] = re_poly->vertex_array[2];
		fourgon->vertex_array[3] = re_poly->vertex_array[3];
		new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if(is_p1_in_tri == 0 && is_p2_in_tri == 0){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[1];
		tri_poly->vertex_array[2] = re_poly->vertex_array[4];
		fourgon->vertex_array[0] = re_poly->vertex_array[1];
		fourgon->vertex_array[1] = re_poly->vertex_array[2];
		fourgon->vertex_array[2] = re_poly->vertex_array[3];
		fourgon->vertex_array[3] = re_poly->vertex_array[4];
		new_triangle_count = trianngulate_4(fourgon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	return new_triangle_count;
}
int trianngulate_6(polygon_t* polygon, int (*points)[2], int** triangles, graph_t* graph, int ind_point,int new_triangle_count){
	polygon_t* re_poly;
	polygon_t* tri_poly;
	polygon_t* fivegon;
	polygon_t* fourgon1;
	polygon_t* fourgon2;
	int orient1, orient2, orient3;
	int i=0;
	int min_x_index = 0;//pi
	int is_p1_in_tri, is_p2_in_tri, is_p3_in_tri;
	min_x_index = find_min_X(polygon,points);
	tri_poly = newPolygon(3);
	fivegon = newPolygon(5);
	fourgon1 = newPolygon(4);
	fourgon2 = newPolygon(4);
	re_poly = newPolygon(polygon->num_vertices);
	//reshuffle points in polygon array
	for(i=0; i<polygon->num_vertices; i++){
		re_poly->vertex_array[i] = polygon->vertex_array[min_x_index];
		min_x_index++;
		if(min_x_index == polygon->num_vertices){
			min_x_index = 0;
		}
	}
	int orient_repoly = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1] );
	if(orient_repoly <0){
		int v = re_poly->vertex_array[1];
		int w = re_poly->vertex_array[2];
		int x = re_poly->vertex_array[3];
		int y = re_poly->vertex_array[4];
		int z= re_poly->vertex_array[5];
		re_poly->vertex_array[1] = z;
		re_poly->vertex_array[2] = y;
		re_poly->vertex_array[3] = x;
		re_poly->vertex_array[4] = w;
		re_poly->vertex_array[5] = v;
	}
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	orient3 = getOrientation(points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[2]][0], points[re_poly->vertex_array[2]][1]);
	is_p1_in_tri = is_point_in_triangle(orient1, orient2, orient3);
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	orient3 = getOrientation(points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[3]][0], points[re_poly->vertex_array[3]][1]);
	is_p2_in_tri = is_point_in_triangle(orient1, orient2, orient3);
	orient1 = getOrientation(points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1]);
	orient2 = getOrientation(points[re_poly->vertex_array[1]][0], points[re_poly->vertex_array[1]][1],points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1]);
	orient3 = getOrientation(points[re_poly->vertex_array[5]][0], points[re_poly->vertex_array[5]][1],points[re_poly->vertex_array[0]][0], points[re_poly->vertex_array[0]][1],points[re_poly->vertex_array[4]][0], points[re_poly->vertex_array[4]][1]);
	is_p3_in_tri = is_point_in_triangle(orient1, orient2, orient3);
	int p1_x, p2_x, p3_x;
	p1_x = points[re_poly->vertex_array[2]][0];
	p2_x = points[re_poly->vertex_array[3]][0];
	p3_x = points[re_poly->vertex_array[4]][0];
	if(is_p1_in_tri==1 && is_p2_in_tri==1 && is_p3_in_tri==1){
		if(p1_x<=p2_x && p1_x<=p3_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[1];
			tri_poly->vertex_array[2] = re_poly->vertex_array[2];
			fivegon->vertex_array[0] = re_poly->vertex_array[0];
			fivegon->vertex_array[1] = re_poly->vertex_array[2];
			fivegon->vertex_array[2] = re_poly->vertex_array[3];
			fivegon->vertex_array[3] = re_poly->vertex_array[4];
			fivegon->vertex_array[4] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_5(fivegon, points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p2_x<=p1_x && p2_x<=p3_x){
			fourgon1->vertex_array[0] = re_poly->vertex_array[0];
			fourgon1->vertex_array[1] = re_poly->vertex_array[1];
			fourgon1->vertex_array[2] = re_poly->vertex_array[2];
			fourgon1->vertex_array[3] = re_poly->vertex_array[3];
			fourgon2->vertex_array[0] = re_poly->vertex_array[0];
			fourgon2->vertex_array[1] = re_poly->vertex_array[3];
			fourgon2->vertex_array[2] = re_poly->vertex_array[4];
			fourgon2->vertex_array[3] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_4(fourgon1,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_4(fourgon2,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p3_x<=p1_x && p3_x<=p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[4];
			tri_poly->vertex_array[2] = re_poly->vertex_array[5];
			fivegon->vertex_array[0] = re_poly->vertex_array[0];
			fivegon->vertex_array[1] = re_poly->vertex_array[1];
			fivegon->vertex_array[2] = re_poly->vertex_array[2];
			fivegon->vertex_array[3] = re_poly->vertex_array[3];
			fivegon->vertex_array[4] = re_poly->vertex_array[4];
			new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
	}
	else if(is_p1_in_tri==0 && is_p2_in_tri==0 && is_p3_in_tri==0){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[1];
		tri_poly->vertex_array[2] = re_poly->vertex_array[5];
		fivegon->vertex_array[0] = re_poly->vertex_array[1];
		fivegon->vertex_array[1] = re_poly->vertex_array[2];
		fivegon->vertex_array[2] = re_poly->vertex_array[3];
		fivegon->vertex_array[3] = re_poly->vertex_array[4];
		fivegon->vertex_array[4] = re_poly->vertex_array[5];
		new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if (is_p1_in_tri == 1 && is_p2_in_tri == 0 && is_p3_in_tri == 0){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[1];
		tri_poly->vertex_array[2] = re_poly->vertex_array[2];
		fivegon->vertex_array[0] = re_poly->vertex_array[0];
		fivegon->vertex_array[1] = re_poly->vertex_array[2];
		fivegon->vertex_array[2] = re_poly->vertex_array[3];
		fivegon->vertex_array[3] = re_poly->vertex_array[4];
		fivegon->vertex_array[4] = re_poly->vertex_array[5];
		new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if (is_p1_in_tri == 0 && is_p2_in_tri == 1 && is_p3_in_tri == 0){
		fourgon1->vertex_array[0] = re_poly->vertex_array[0];
		fourgon1->vertex_array[1] = re_poly->vertex_array[1];
		fourgon1->vertex_array[2] = re_poly->vertex_array[2];
		fourgon1->vertex_array[3] = re_poly->vertex_array[3];
		fourgon2->vertex_array[0] = re_poly->vertex_array[0];
		fourgon2->vertex_array[1] = re_poly->vertex_array[3];
		fourgon2->vertex_array[2] = re_poly->vertex_array[4];
		fourgon2->vertex_array[3] = re_poly->vertex_array[5];
		new_triangle_count = trianngulate_4(fourgon1,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_4(fourgon2,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if (is_p1_in_tri == 0 && is_p2_in_tri == 0 && is_p3_in_tri == 1){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[4];
		tri_poly->vertex_array[2] = re_poly->vertex_array[5];
		fivegon->vertex_array[0] = re_poly->vertex_array[0];
		fivegon->vertex_array[1] = re_poly->vertex_array[1];
		fivegon->vertex_array[2] = re_poly->vertex_array[2];
		fivegon->vertex_array[3] = re_poly->vertex_array[3];
		fivegon->vertex_array[4] = re_poly->vertex_array[4];
		new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if (is_p1_in_tri == 1 && is_p2_in_tri == 1 && is_p3_in_tri == 0){
		if(p1_x<p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[1];
			tri_poly->vertex_array[2] = re_poly->vertex_array[2];
			fivegon->vertex_array[0] = re_poly->vertex_array[0];
			fivegon->vertex_array[1] = re_poly->vertex_array[2];
			fivegon->vertex_array[2] = re_poly->vertex_array[3];
			fivegon->vertex_array[3] = re_poly->vertex_array[4];
			fivegon->vertex_array[4] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p1_x>p2_x){
			fourgon1->vertex_array[0] = re_poly->vertex_array[0];
			fourgon1->vertex_array[1] = re_poly->vertex_array[1];
			fourgon1->vertex_array[2] = re_poly->vertex_array[2];
			fourgon1->vertex_array[3] = re_poly->vertex_array[3];
			fourgon2->vertex_array[0] = re_poly->vertex_array[0];
			fourgon2->vertex_array[1] = re_poly->vertex_array[3];
			fourgon2->vertex_array[2] = re_poly->vertex_array[4];
			fourgon2->vertex_array[3] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_4(fourgon1,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_4(fourgon2,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p1_x == p2_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[1];
			tri_poly->vertex_array[2] = re_poly->vertex_array[2];
			fivegon->vertex_array[0] = re_poly->vertex_array[0];
			fivegon->vertex_array[1] = re_poly->vertex_array[2];
			fivegon->vertex_array[2] = re_poly->vertex_array[3];
			fivegon->vertex_array[3] = re_poly->vertex_array[4];
			fivegon->vertex_array[4] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
	}
	else if (is_p1_in_tri == 1 && is_p2_in_tri == 0 && is_p3_in_tri == 1){
		tri_poly->vertex_array[0] = re_poly->vertex_array[0];
		tri_poly->vertex_array[1] = re_poly->vertex_array[1];
		tri_poly->vertex_array[2] = re_poly->vertex_array[2];
		fivegon->vertex_array[0] = re_poly->vertex_array[0];
		fivegon->vertex_array[1] = re_poly->vertex_array[2];
		fivegon->vertex_array[2] = re_poly->vertex_array[3];
		fivegon->vertex_array[3] = re_poly->vertex_array[4];
		fivegon->vertex_array[4] = re_poly->vertex_array[5];
		new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
		new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
	}
	else if (is_p1_in_tri == 0 && is_p2_in_tri == 1 && is_p3_in_tri == 1){
		if(p2_x<p3_x){
			fourgon1->vertex_array[0] = re_poly->vertex_array[0];
			fourgon1->vertex_array[1] = re_poly->vertex_array[1];
			fourgon1->vertex_array[2] = re_poly->vertex_array[2];
			fourgon1->vertex_array[3] = re_poly->vertex_array[3];
			fourgon2->vertex_array[0] = re_poly->vertex_array[0];
			fourgon2->vertex_array[1] = re_poly->vertex_array[3];
			fourgon2->vertex_array[2] = re_poly->vertex_array[4];
			fourgon2->vertex_array[3] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_4(fourgon1,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_4(fourgon2,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p2_x > p3_x){
			tri_poly->vertex_array[0] = re_poly->vertex_array[0];
			tri_poly->vertex_array[1] = re_poly->vertex_array[4];
			tri_poly->vertex_array[2] = re_poly->vertex_array[5];
			fivegon->vertex_array[0] = re_poly->vertex_array[0];
			fivegon->vertex_array[1] = re_poly->vertex_array[1];
			fivegon->vertex_array[2] = re_poly->vertex_array[2];
			fivegon->vertex_array[3] = re_poly->vertex_array[3];
			fivegon->vertex_array[4] = re_poly->vertex_array[4];
			new_triangle_count = trianngulate_5(fivegon,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_3(tri_poly,points,triangles,graph,ind_point,new_triangle_count);
		}
		else if(p2_x == p3_x){
			fourgon1->vertex_array[0] = re_poly->vertex_array[0];
			fourgon1->vertex_array[1] = re_poly->vertex_array[1];
			fourgon1->vertex_array[2] = re_poly->vertex_array[2];
			fourgon1->vertex_array[3] = re_poly->vertex_array[3];
			fourgon2->vertex_array[0] = re_poly->vertex_array[0];
			fourgon2->vertex_array[1] = re_poly->vertex_array[3];
			fourgon2->vertex_array[2] = re_poly->vertex_array[4];
			fourgon2->vertex_array[3] = re_poly->vertex_array[5];
			new_triangle_count = trianngulate_4(fourgon1,points,triangles,graph,ind_point,new_triangle_count);
			new_triangle_count = trianngulate_4(fourgon2,points,triangles,graph,ind_point,new_triangle_count);
		}
	}
	return new_triangle_count;
}

ploc_t* create_ploc(int (*points)[2], int (*triangles)[3], int n, int m){
	int *indep_set;
	int *indep_set_points;
    int *ind_set_touch_tri_flag;
	int i, j, v;
	int a,b,c;
	int **triangle_new;
	int seen_count;
	adj_list_node_t *pCrawl;
	int ind_set_count = 0;
    polygon_t* poly;
	int triangle_index;
	int ind_point;
	int prev;
	int new_triangle_count = 0;
	graph_t* graph = NULL;
	ploc_t* root = NULL;
	ploc_t* tail = NULL;
	ploc_t* pl = NULL;
	ploc_t* tmp = NULL;
	pl = (ploc_t*)malloc(sizeof(ploc_t));
	pl->points = (int**) malloc(n * sizeof(int*));
	if((pl->points) == NULL){
		fprintf(stderr, "out of memory\n");
		exit(0);
	}
	for(i = 0; i < n; i++){
		(pl->points)[i] = (int*) malloc(2 * sizeof(int));
		if((pl->points)[i] == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
	}
	for(i=0; i<n;i++){
		for(j=0; j<2;j++){
			(pl->points)[i][j] = points[i][j];
		}
	}
	pl->triangle = (int**) malloc(m * sizeof(int*));
	if((pl->triangle) == NULL){
		fprintf(stderr, "out of memory\n");
		exit(0);
	}
	for(i = 0; i < m; i++){
		(pl->triangle)[i] = (int*) malloc(3 * sizeof(int));
		if((pl->triangle)[i] == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
	}
	for(i=0; i<m;i++){
		for(j=0; j<3;j++){
			(pl->triangle)[i][j] = triangles[i][j];
		}
	}
	pl->previous = NULL;
	pl->next = NULL;
	pl->tri_count = m;
	root  = pl;
	tail=pl;
	while(tail->tri_count>=3){
		ind_set_count = 0;
		new_triangle_count = 0;
		m = tail->tri_count;
		triangle_new = (int**) malloc(m * sizeof(int*));
		if(triangle_new == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
		for(i = 0; i < m; i++){
			triangle_new[i] = (int*) malloc(9 * sizeof(int));
			if(triangle_new[i] == NULL){
				fprintf(stderr, "out of memory\n");
				exit(0);
			}
		}
		for(i = 0; i < m; i++){
			for(j = 0; j < 9; j++){
				triangle_new[i][j] = -1;
			}
		}
		indep_set = (int*) malloc((n) * sizeof(int));
		if(indep_set == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
		for(i = 0; i < n; i++){
			indep_set[i] = -1;
		}
		ind_set_touch_tri_flag = (int*)calloc(m, sizeof(int));
		if(ind_set_touch_tri_flag == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
		graph = createGraph(n);
		for(i=0; i<m;i++){
			a = (tail->triangle)[i][0];
			b = (tail->triangle)[i][1];
			c = (tail->triangle)[i][2];
			(graph->array[a]).pt_degree = (((graph->array[a]).pt_degree) + 1);
			(graph->array[b]).pt_degree = (((graph->array[b]).pt_degree) + 1);
			(graph->array[c]).pt_degree = (((graph->array[c]).pt_degree) + 1);
		}
		for(i=0; i<m;i++){
			a = (tail->triangle)[i][0];
			b = (tail->triangle)[i][1];
			c = (tail->triangle)[i][2];
			addPointTriangle(graph,a,i);
			addPointTriangle(graph,b,i);
			addPointTriangle(graph,c,i);
		}
		indep_set[90000] = indep_set[90001] = indep_set[90002] = 0;
		for(i=0; i<n;i++){
			if(((graph->array[i]).pt_degree) > 6){
				indep_set[i] = 0;
				graph->array[i].in_indep_set = 0;
			}
			else if(indep_set[i] != 0 &&(((graph->array[i]).pt_degree) !=0)){
				indep_set[i] = 1;
				graph->array[i].in_indep_set = 1;
				ind_set_count++;
				pCrawl = graph->array[i].head;
				while(pCrawl){
					a = (tail->triangle)[pCrawl->triangle_no][0];
					b = (tail->triangle)[pCrawl->triangle_no][1];
					c = (tail->triangle)[pCrawl->triangle_no][2];
					graph->array[a].in_indep_set = 0;
					graph->array[b].in_indep_set = 0;
					graph->array[c].in_indep_set = 0;
					indep_set[a] = 0;
					indep_set[b] = 0;
					indep_set[c] = 0;
					ind_set_touch_tri_flag[pCrawl->triangle_no] = 1;
					pCrawl = pCrawl->next;
				}
				indep_set[i] = 1;
				graph->array[i].in_indep_set = 1;
			}
		}
		//printf("ind_set_count=%d\n",ind_set_count);
		indep_set_points = (int*)malloc((ind_set_count) * sizeof(int));
		if(indep_set_points == NULL){
				fprintf(stderr, "out of memory\n");
				exit(0);
		}
		j=0;
		for(i = 0; i < n; i++){
			if(indep_set[i] == 1)
				indep_set_points[j++] = i;
		}
		for(i=0; i<m;i++){
			if(ind_set_touch_tri_flag[i] == 0){
				a = (tail->triangle)[i][0];
				b = (tail->triangle)[i][1];
				c = (tail->triangle)[i][2];
				triangle_new[new_triangle_count][0] = a;
				triangle_new[new_triangle_count][1] = b;
				triangle_new[new_triangle_count][2] = c;
				triangle_new[new_triangle_count][3] = i;
				triangle_new[new_triangle_count][4] = -1;
				triangle_new[new_triangle_count][5] = -1;
				triangle_new[new_triangle_count][6] = -1;
				triangle_new[new_triangle_count][7] = -1;
				triangle_new[new_triangle_count][8] = -1;
				new_triangle_count++;
			}
		}
		for(i=0; i<ind_set_count;i++){
			ind_point = indep_set_points[i];
			pCrawl = graph->array[ind_point].head;
			v = graph->array[ind_point].pt_degree;
			poly = newPolygon(v);
			seen_count = 0;
			j=0;
			triangle_index = pCrawl->triangle_no;
			a = (tail->triangle)[triangle_index][0];
			b = (tail->triangle)[triangle_index][1];
			c = (tail->triangle)[triangle_index][2];
			if(a == ind_point){
				poly->vertex_array[j++] = b;
				poly->vertex_array[j++] = c;
				prev = c;
			}
			else if(b == ind_point){
				poly->vertex_array[j++] = a;
				poly->vertex_array[j++] = c;
				prev = c;
			}
			else{
				poly->vertex_array[j++] = a;
				poly->vertex_array[j++] = b;
				prev = b;
			}
			pCrawl->visited = 1;
			seen_count++;
			pCrawl = pCrawl->next;
			while(seen_count!=v && j!=v){
				if(pCrawl == NULL)
					pCrawl = graph->array[ind_point].head;
				if(pCrawl->visited ==0){
					triangle_index = pCrawl->triangle_no;
					a = (tail->triangle)[triangle_index][0];
					b = (tail->triangle)[triangle_index][1];
					c = (tail->triangle)[triangle_index][2];
					if(prev == a || prev == b || prev == c){
						seen_count++;
						pCrawl->visited = 1;
						if(a == ind_point){
							if(prev==b){
								poly->vertex_array[j++] = c;
								prev = c;
							}
							else if(prev==c){
								poly->vertex_array[j++] = b;
								prev = b;
							}
						}
						else if(b == ind_point){
							if(prev==a){
								poly->vertex_array[j++] = c;
								prev = c;
							}
							else{
								poly->vertex_array[j++] = a;
								prev = a;
							}
						}
						else{
							if(prev==a){
								poly->vertex_array[j++] = b;
								prev = b;
							}
							else{
								poly->vertex_array[j++] = a;
								prev = a;
							}
						}
						//printf("j=%d\n",j);
					}
				}

				pCrawl = pCrawl->next;
			}

			//retriangulate call functions based on number of vertices.
			if(poly->num_vertices == 3){
				new_triangle_count = trianngulate_3(poly, points,triangle_new, graph, ind_point, new_triangle_count);
			}
			else if(poly->num_vertices == 4){
				new_triangle_count = trianngulate_4(poly, points, triangle_new, graph, ind_point, new_triangle_count);
			}
			else if(poly->num_vertices == 5){
				new_triangle_count = trianngulate_5(poly,points, triangle_new, graph, ind_point, new_triangle_count);
			}
			else if(poly->num_vertices == 6){
				new_triangle_count = trianngulate_6(poly,points, triangle_new, graph, ind_point, new_triangle_count);
			}
		}
		tmp = (ploc_t*)malloc(sizeof(ploc_t));
		tmp->points = (int**) malloc(n * sizeof(int*));
		if((tmp->points) == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
		for(i = 0; i < n; i++){
			(tmp->points)[i] = (int*) malloc(2 * sizeof(int));
			if((tmp->points)[i] == NULL){
				fprintf(stderr, "out of memory\n");
				exit(0);
			}
		}
		for(i=0; i<n;i++){
			for(j=0; j<2;j++){
				(tmp->points)[i][j] = points[i][j];
			}
		}
		tmp->triangle = (int**) malloc(new_triangle_count * sizeof(int*));
		if((tmp->triangle) == NULL){
			fprintf(stderr, "out of memory\n");
			exit(0);
		}
		for(i = 0; i < new_triangle_count; i++){
			(tmp->triangle)[i] = (int*) malloc(9 * sizeof(int));
			if((tmp->triangle)[i] == NULL){
				fprintf(stderr, "out of memory\n");
				exit(0);
			}
		}
		for(i=0; i<new_triangle_count;i++){
			for(j=0; j<9;j++){
				(tmp->triangle)[i][j] = triangle_new[i][j];
			}
		}
		tmp->previous = tail;
		tail = tmp;
		tmp->previous->next = tmp;
		tmp->next = NULL;
		tmp->tri_count = new_triangle_count;
	//}
	}
	return tail;
}
int query_ploc(ploc_t* pl,int x, int y){
	ploc_t* tmp = pl;
	int i,j,k,orient1,orient2,orient3,is_pt_in_tri;
	int** t = tmp->triangle;
	int** p = tmp->points;
	int* tp = NULL;
	int sub_triangle = -1;
	orient1 = getOrientation(p[t[1][0]][0],p[t[1][0]][1],p[t[1][1]][0],p[t[1][1]][1],x,y);
	orient2 = getOrientation(p[t[1][1]][0],p[t[1][1]][1],p[t[1][2]][0],p[t[1][2]][1],x,y);
	orient3 = getOrientation(p[t[1][2]][0],p[t[1][2]][1],p[t[1][0]][0],p[t[1][0]][1],x,y);
	is_pt_in_tri = is_point_in_triangle(orient1,orient2, orient3);
	if(is_pt_in_tri == 0){
		return -1;
	}
	else{
		tp = (int*) malloc(6*sizeof(int));
		for(i=0; i<6;i++){
			tp[i] = -1;
		}
		j=0;
		for(i=3;i<9;i++){
			tp[j++] = t[1][i];
		}
		tmp = tmp->previous;
		t = tmp->triangle;
		p = tmp->points;
		while(tmp->previous != NULL){
			for(i=0;i<6;i++){
				if(tp[i]!=-1){
					orient1 = getOrientation(p[t[tp[i]][0]][0],p[t[tp[i]][0]][1],p[t[tp[i]][1]][0],p[t[tp[i]][1]][1],x,y);
					orient2 = getOrientation(p[t[tp[i]][1]][0],p[t[tp[i]][1]][1],p[t[tp[i]][2]][0],p[t[tp[i]][2]][1],x,y);
					orient3 = getOrientation(p[t[tp[i]][2]][0],p[t[tp[i]][2]][1],p[t[tp[i]][0]][0],p[t[tp[i]][0]][1],x,y);
					is_pt_in_tri = is_point_in_triangle(orient1,orient2, orient3);
					if(is_pt_in_tri == 0){
						continue;
					}
					else{
						sub_triangle = tp[i];
						for(j=0; j<6;j++){
							tp[j] = -1;
						}
						j=0;
						for(k=3;k<9;k++){
							tp[j++] = t[sub_triangle][k];
						}
						tmp = tmp->previous;
						t = tmp->triangle;
						p = tmp->points;
						break;
					}
				}
			}
		}
		for(i=0;i<6;i++){
			if(tp[i]!=-1){
				orient1 = getOrientation(p[t[tp[i]][0]][0],p[t[tp[i]][0]][1],p[t[tp[i]][1]][0],p[t[tp[i]][1]][1],x,y);
				orient2 = getOrientation(p[t[tp[i]][1]][0],p[t[tp[i]][1]][1],p[t[tp[i]][2]][0],p[t[tp[i]][2]][1],x,y);
				orient3 = getOrientation(p[t[tp[i]][2]][0],p[t[tp[i]][2]][1],p[t[tp[i]][0]][0],p[t[tp[i]][0]][1],x,y);
				is_pt_in_tri = is_point_in_triangle(orient1,orient2, orient3);
				if(is_pt_in_tri == 0){
					sub_triangle = -1;
					continue;
				}
				else{
					sub_triangle = tp[i];
					for(j=0; j<6;j++){
						tp[j] = -1;
					}
					j=0;
					for(k=3;k<9;k++){
						tp[j++] = t[sub_triangle][k];
					}
					tmp = tmp->previous;
					break;
				}
			}
		}
		return sub_triangle;

	}
}
int main() {
	int points[90003][2];
	int triangles[180002][3];
	int i, j, k;	k = 0;
	for (i = 0; i < 300; i++) {
		for (j = 0; j < 300; j++) {
			points[k][0] = 15 * i;
			points[k][1] = 15 * j;
			//printf("points[%d][x]=%d, points[%d][y]=%d\n",k,points[k][0],k,points[k][1]);
			k += 1;
		}
	}
	points[k][0] = -4500;
	points[k++][1] = -15;
	points[k][0] = 9000;
	points[k++][1] = -15;
	points[k][0] = 2250;
	points[k++][1] = 9000;
	printf("Prepared %d points\n", k);
	k = 0;
	triangles[k][0] = 90000;
	triangles[k][1] = 90001;
	triangles[k++][2] = 90002;
	for (i = 0; i < 299; i++) {
		triangles[k][0] = i;
		triangles[k][1] = i + 1;
		triangles[k++][2] = 90000;
	}
	triangles[k][0] = 299;
	triangles[k][1] = 90000;
	triangles[k++][2] = 90002;
	for (i = 0; i < 299; i++) {
		triangles[k][0] = 299 + 300 * i;
		triangles[k][1] = 299 + 300 * (i + 1);
		triangles[k++][2] = 90002;
	}
	triangles[k][0] = 90002;
	triangles[k][1] = 89999;
	triangles[k++][2] = 90001;
	for (i = 0; i < 299; i++) {
		triangles[k][0] = 299 * 300 + i;
		triangles[k][1] = 299 * 300 + i + 1;
		triangles[k++][2] = 90001;
	}
	triangles[k][0] = 90000;
	triangles[k][1] = 89700;
	triangles[k++][2] = 90001;
	for (i = 0; i < 299; i++) {
		triangles[k][0] = 300 * i;
		triangles[k][1] = 300 * (i + 1);
		triangles[k++][2] = 90000;
	}
	for (i = 0; i < 299; i++) {
		for (j = 0; j < 299; j++) {
			triangles[k][0] = 300 * i + j;
			triangles[k][1] = 300 * (i + 1) + j;
			triangles[k++][2] = 300 * i + j + 1;
			triangles[k][0] = 300 * (i + 1) + j + 1;
			triangles[k][1] = 300 * (i + 1) + j;
			triangles[k++][2] = 300 * i + j + 1;
		}
	}
	printf("Prepared %d triangles\n", k);

	ploc_t *pl;
	pl = create_ploc(points, triangles, 90003, 180002);

	for (i = 0; i < 10000; i++) {
		int a, b, c, x, y, t;
		j = (rand() % 180001) + 1;
		a = triangles[j][0];
		b = triangles[j][1];
		c = triangles[j][2];
		x = (points[a][0] + points[b][0] + points[c][0]) / 3;
		y = (points[a][1] + points[b][1] + points[c][1]) / 3;
		t = query_ploc(pl, x, y);
		if (t != j) {
			printf("Error on triangle %d, misidentified as triangle %d\n",
					j, t);
			printf(
					"Point (%d,%d) should be in triangle (%d,%d), (%d,%d), (%d,%d)\n",
					x, y, points[a][0], points[a][1], points[b][0],
					points[b][1], points[c][0], points[c][1]);
			printf(
					"Instead claimed in triangle (%d,%d), (%d,%d), (%d,%d)\n",
					points[triangles[t][0]][0], points[triangles[t][0]][1],
					points[triangles[t][1]][0], points[triangles[t][1]][1],
					points[triangles[t][2]][0], points[triangles[t][2]][1]);
			exit(0);
		}
	}

	printf("Passed test\n");
}
