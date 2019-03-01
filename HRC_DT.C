#include </home/kristi/Documents/research/cleap/src/cleap.h>
#include </home/kristi/Documents/research/cleap/src/cleap_private.h>

#include <map>
#include <string>
#include <tr1/unordered_map>
#include <iostream>
#include "visit_writer.c"
#include <vector>

using std::vector;

using std::cerr;
using std::endl;

#define NUMPOINTS 10
#define DIM 2

struct Pair {
  float x;
  float y;
};

//Generates random data point values
void PointGenerator(Pair* DTarray) {
  for(int j = 0; j < NUMPOINTS; j++) {

    float rand_value = rand() % 1000000 / 1000000.0;
    DTarray[j].x = rand_value;

    rand_value = rand() % 1000000 / 1000000.0;
    DTarray[j].y = rand_value;
  }

}

//compare function for qsort - sorts along X axis first
int compareX (const void * a, const void * b)
{
  const Pair* A = (const Pair*) a;
  const Pair* B = (const Pair*) b;
  if (A->x > B->x) return 1;
  else if (A->x < B->x) return -1;
  else if (A->x == B->x) {
    if (A->y > B->y) return 1;
    else if (A->y < B->y) return -1;
    else if (A->y == B->y) return 0;
  }
  else return EXIT_FAILURE;

  return 0;
}

//Partition the data points into slabs - function defines one slab
void slabPartition(Pair *DTarray, float *slabs, int factor, int offset, const int pop) {
  int count = 0; int j = 0;
  while (count != pop) {
    slabs[j] = DTarray[offset + count].x;
    slabs[j + 1] = DTarray[offset + count].y;
    j+=2; count++;
  }
}

bool IsOnSameSide(float *endPoint1, float *endPoint2, 
                  float *referencePoint, float *newPoint)
{

    // see: http://doubleroot.in/lessons/straight-line/position-of-a-point-relative-to-a-line/#.Wt5H7ZPwalM


    float m, b;
    // need to solve equation y = mx + b for endPoint1 
    // and endPoint2.

    if (endPoint1[0] == endPoint2[0])
    {
        // infinite slope ... fail
        return false;
    }
    m = (endPoint2[1] - endPoint1[1])/(endPoint2[0] - endPoint1[0]);
    // y = mx+b
    // a'x+b'y+c' = 0
    // mx-y+b = 0;
    // a' = m, b' = -1, c' = b
    b = endPoint2[1]-m*endPoint2[0];
    float a_formula = m;
    float b_formula = -1;
    float c_formula = b;

    float val1 = referencePoint[0]*a_formula + referencePoint[1]*b_formula + c_formula;
    float val2 = newPoint[0]*a_formula + newPoint[1]*b_formula + c_formula;

    float product = val1*val2;
    return (product < 0 ? false : true);
}

class OneTriangle
{
  public:
    float     p1[2]; 
    float     p2[2]; 
    float     p3[2]; 

    bool      ContainsPoint(float x, float y);
};

bool
OneTriangle::ContainsPoint(float x, float y)
{
    float p4[2];
    p4[0] = x;
    p4[1] = y;
    bool p3_and_p4 = IsOnSameSide(p1, p2, p3, p4);
    bool p1_and_p4 = IsOnSameSide(p3, p2, p1, p4);
    bool p2_and_p4 = IsOnSameSide(p3, p1, p2, p4);
    if (p3_and_p4 && p1_and_p4 && p2_and_p4)
        return true;
    return false;
}

class DelaunayTriangulation
{
  public:
    DelaunayTriangulation();
    ~DelaunayTriangulation();
    void   Initialize(float, float, float, float, float, float);
    void   AddPoint(float, float);
    int    TriGetSize() {return triangles.size();};
    int    FindFaces(float*, const int);
    void   VerifyResults(float*, const int);    
    void   Clear();
    //void   WriteOutTriangle(char *filename);

    int                       count;
    std::vector<OneTriangle>  triangles;
    std::vector<int>          faces;
};
/*
void DelaunayTriangulation::WriteOutTriangle(char *filename)
{
    int ncells = triangles.size();
cerr << "NUMBER OF TRIANGLE is " << ncells << endl;
    int *celltypes = new int[ncells];
    for (int i = 0 ; i < ncells ; i++)
        celltypes[i] = VISIT_TRIANGLE;

    int dimensions = 3; // always 3 for VTK
    int vertices_per_cell = 3;
    int npts = ncells*vertices_per_cell*dimensions;
    float *pts = new float[npts];
    int *conn = new int[ncells*vertices_per_cell];
    int offset = 0;
    for (int i = 0 ; i < ncells ; i++)
    {
        pts[offset+0] = triangles[i].p1[0];
        pts[offset+1] = triangles[i].p1[1];
        pts[offset+2] = 0;
        offset += 3;
        pts[offset+0] = triangles[i].p2[0];
        pts[offset+1] = triangles[i].p2[1];
        pts[offset+2] = 0;
        offset += 3;
        pts[offset+0] = triangles[i].p3[0];
        pts[offset+1] = triangles[i].p3[1];
        pts[offset+2] = 0;
        offset += 3;
    }

    for (int i = 0 ; i < 3*ncells ; i++)
    {
        conn[i] = i;
    }
    write_unstructured_mesh(filename, 0, npts/3, pts,
                            ncells, celltypes, conn, 0,
                            NULL, NULL, NULL, NULL);
}  
*/

//constructor to initialize count private variable
DelaunayTriangulation::DelaunayTriangulation() { count = 0; }

//destructor to clear out memory
DelaunayTriangulation::~DelaunayTriangulation() {
  Clear();
  triangles.shrink_to_fit(); 
  faces.shrink_to_fit();
}

void
DelaunayTriangulation::Initialize(float x1, float y1, float x2, float y2, float x3, float y3)
{
    OneTriangle ot;
    ot.p1[0] = x1;
    ot.p1[1] = y1;
    ot.p2[0] = x2;
    ot.p2[1] = y2;
    ot.p3[0] = x3;
    ot.p3[1] = y3;
    triangles.push_back(ot);
}

void
DelaunayTriangulation::AddPoint(float x1, float y1)
{
  for (int i = 0 ; i < triangles.size() ; i++)
  {
    if (triangles[i].ContainsPoint(x1, y1))
    {
      OneTriangle original_triangle = triangles[i];
      // split triangle i into three triangles
      // note: no edge flipping or Delaunay business.
      // start by replacing triangle in the current list 
      triangles[i].p3[0] = x1;
      triangles[i].p3[1] = y1; 

      // now add two more triangles.
      OneTriangle new_triangle1;
      new_triangle1.p1[0] = x1;
      new_triangle1.p1[1] = y1;
      new_triangle1.p2[0] = original_triangle.p2[0];
      new_triangle1.p2[1] = original_triangle.p2[1]; 
      new_triangle1.p3[0] = original_triangle.p3[0];
      new_triangle1.p3[1] = original_triangle.p3[1];
      triangles.push_back(new_triangle1);

      OneTriangle new_triangle2;
      new_triangle2.p1[0] = original_triangle.p1[0];
      new_triangle2.p1[1] = original_triangle.p1[1];
      new_triangle2.p2[0] = x1;
      new_triangle2.p2[1] = y1;
      new_triangle2.p3[0] = original_triangle.p3[0];
      new_triangle2.p3[1] = original_triangle.p3[1];
      triangles.push_back(new_triangle2);

      break;
    }
  }
}

int
DelaunayTriangulation::FindFaces(float* slabs, const int pop) { 
  int size = triangles.size();

  for(int i = 0; i < size; i++) {   
    for(int s = 0; s < (pop+3); s++) {
      if(slabs[2*s] == triangles[i].p1[0] && slabs[2*s +1] == triangles[i].p1[1]) {
        faces.push_back(s);
      }
      if(slabs[2*s] == triangles[i].p2[0] && slabs[2*s +1] == triangles[i].p2[1]) {
	faces.push_back(s);
      }
      if(slabs[2*s] == triangles[i].p3[0] && slabs[2*s +1] == triangles[i].p3[1]) {
        faces.push_back(s);	
      }
    } 
  } 
}

//to double check results
void 
DelaunayTriangulation::VerifyResults (float* slabs, const int pop) {
  for(int i = 0; i < triangles.size(); i++) {
    printf("Triangle %d: %lf, %lf\n %lf, %lf\n %lf, %lf\n", i+1, triangles[i].p1[0], triangles[i].p1[1], triangles[i].p2[0], triangles[i].p2[1], triangles[i].p3[0], triangles[i].p3[1]);
  }

  for(int i = 0; i < (pop + 3); i++) {
    printf("\nFor vertex %d: xcoord: %lf, ycoord: %lf \n", i, slabs[2*i], slabs[2*i+1]);
  }

  for(int i = 0; i < faces.size(); i+=3) {
    printf("face %d, %d, %d\n", faces[i], faces[i+1], faces[i+2]);
  }

  printf("\n");
}

void
DelaunayTriangulation::Clear() {
  if(triangles.size() > 0)  
    triangles.clear();
  if(faces.size() > 0)
    faces.clear();
}
/*************************BEGIN CLEAP FUNCTIONS*******************/

CLEAP_RESULT _cleap_generate_edges_hash(_cleap_mesh *m, DelaunayTriangulation dt){
  //parsing faces and edges
  int face = 3;
  float3 normal;
  float3 v1,v2;

  // the hash for edges
  std::tr1::unordered_map<int, std::tr1::unordered_map<int, _tmp_edge> > root_hash;
  std::tr1::unordered_map<int, _tmp_edge>::iterator hit;
  std::vector<_tmp_edge*> edge_vector;

  _tmp_edge* aux_tmp_edge;
  int j_sec[3] = {0, 0, 1};
  int k_sec[3] = {1, 2, 2};
  int op_sec[3] = {2, 1, 0};
  int j ,k, op;
 
  int count = 0; //do i need this?

  for(int i=0; i<m->face_count; i++) {  
    m->triangles[i*3] = dt.faces[i*3];
    m->triangles[i*3+1] = dt.faces[i*3+1];
    m->triangles[i*3+2] = dt.faces[i*3+2];

    //count+=3;

    //Building Edges
    for(int q=0; q<3; q++){
      j=j_sec[q], k=k_sec[q], op=op_sec[q];
      // always the higher first
      if( m->triangles[i*3+j] < m->triangles[i*3+k]){
        k = j;
        j = k_sec[q];
      }
    }

    // ok, first index already existed, check if the second exists or not
    std::tr1::unordered_map<int, _tmp_edge> *second_hash = &root_hash[m->triangles[i*3+j]];
    hit = second_hash->find(m->triangles[i*3+k]);
    if( hit != second_hash->end() ){
      // the edge already exists, then fill the remaining info
      aux_tmp_edge = &(hit->second);
      aux_tmp_edge->b1 = i*3+j;
      aux_tmp_edge->b2 = i*3+k;
      aux_tmp_edge->op2 = i*3+op;
    } else{
      // create a new edge
      aux_tmp_edge = &(*second_hash)[m->triangles[i*3+k]]; // create the low value on secondary_hash
      aux_tmp_edge->n1 = m->triangles[i*3+j];
      aux_tmp_edge->n2 = m->triangles[i*3+k];
      aux_tmp_edge->a1 = i*3+j;
      aux_tmp_edge->a2 = i*3+k;
      aux_tmp_edge->b1 = -1;
      aux_tmp_edge->b2 = -1;
      aux_tmp_edge->op1 = i*3+op;
      aux_tmp_edge->op2 = -1;

      aux_tmp_edge->id = edge_vector.size();
      edge_vector.push_back( aux_tmp_edge );
    }

    float4 p1 = m->vnc_data.v[m->triangles[i*face]];
    float4 p2 = m->vnc_data.v[m->triangles[i*face+1]];
    float4 p3 = m->vnc_data.v[m->triangles[i*face+2]];
    v1 = make_float3( p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
    v2 = make_float3( p3.x - p1.x, p3.y - p1.y, p3.z - p1.z);
    normal.x =   (v1.y * v2.z) - (v2.y * v1.z);
    normal.y = -((v1.x * v2.z) - (v2.x * v1.z));
    normal.z =   (v1.x * v2.y) - (v2.x * v1.y);
    //printf("   Normal= (%f, %f, %f)\n", normal.x, normal.y, normal.z);

    //!Calculate Normals for this face
    m->vnc_data.n[m->triangles[i*face]].x += normal.x;
    m->vnc_data.n[m->triangles[i*face]].y += normal.y;
    m->vnc_data.n[m->triangles[i*face]].z += normal.z;

    m->vnc_data.n[m->triangles[i*face+1]].x += normal.x;
    m->vnc_data.n[m->triangles[i*face+1]].y += normal.y;
    m->vnc_data.n[m->triangles[i*face+1]].z += normal.z;

    m->vnc_data.n[m->triangles[i*face+2]].x += normal.x;
    m->vnc_data.n[m->triangles[i*face+2]].y += normal.y;
    m->vnc_data.n[m->triangles[i*face+2]].z += normal.z;
  }

  m->processed_edges = 0;
  // CLEAP::MESH:: update the edge count, now after being calculated
  m->edge_count = edge_vector.size();
  // CLEAP::MESH:: malloc edge data
  m->edge_data.n = (int2*)malloc( sizeof(int2)*m->edge_count );
  m->edge_data.a = (int2*)malloc( sizeof(int2)*m->edge_count );
  m->edge_data.b = (int2*)malloc( sizeof(int2)*m->edge_count );
  m->edge_data.op = (int2*)malloc( sizeof(int2)*m->edge_count );
  // CLEAP::MESH:: put edge data into its final format that matches the _cleap_mesh structure

  for( int i=0; i<m->edge_count; i++ ){
    m->edge_data.n[i] = make_int2(edge_vector[i][0].n1, edge_vector[i][0].n2);
    m->edge_data.a[i] = make_int2(edge_vector[i][0].a1, edge_vector[i][0].a2);
    m->edge_data.b[i] = make_int2(edge_vector[i][0].b1, edge_vector[i][0].b2);
    m->edge_data.op[i] = make_int2(edge_vector[i][0].op1, edge_vector[i][0].op2);
  }

  edge_vector.clear();
  //printf("ok\n"); fflush(stdout);
}

CLEAP_RESULT load_mesh_host(const int pop, float* slabs, _cleap_mesh *m, DelaunayTriangulation dt) {
  //prepare cleap mesh
  m->vertex_count = pop+3; //add 3 to include boundary triangle vertices
  m->face_count = dt.TriGetSize(); //faces is equal to number of triangles in mesh
  m->edge_count = (pop+3) + (dt.TriGetSize()) - 2; //following Euler formula, V-E+F=2, so E=V+F-2

  _cleap_reset_minmax(m); //do i need this?

  m->triangles = (GLuint*)malloc(sizeof(GLuint)*m->face_count*3); //malloc host triangles array
  m->vnc_data.v = (float4*)malloc(sizeof(float4)*m->vertex_count); //malloc vertex data => struct of arrays
  m->vnc_data.n = (float4*)malloc(sizeof(float4)*m->vertex_count); //malloc vertex data => struct of arrays
  m->vnc_data.c = (float4*)malloc(sizeof(float4)*m->vertex_count); //malloc vertex data => struct of arrays

  //parse vertex data
  int count = 0;
  for(int i = 0; i < m->vertex_count; i++) {
    m->vnc_data.v[count].x = slabs[2*i];
    m->vnc_data.v[count].y = slabs[2*i+1];
    m->vnc_data.v[count].z = 0; //only doing 2D for now
    count++;
   
    m->vnc_data.v[i].w = 1.0f;

    // normals
    m->vnc_data.n[i] = make_float4(0.0f, 0.0f, 0.0f, 1.0f);

    // maximum values
    if (m->vnc_data.v[i].x > m->max_coords.x) m->max_coords.x=m->vnc_data.v[i].x;
    if (m->vnc_data.v[i].y > m->max_coords.y) m->max_coords.y=m->vnc_data.v[i].y;
    if (m->vnc_data.v[i].z > m->max_coords.z) m->max_coords.z=m->vnc_data.v[i].z;
    if (m->vnc_data.v[i].x < m->min_coords.x) m->min_coords.x=m->vnc_data.v[i].x;
    if (m->vnc_data.v[i].y < m->min_coords.y) m->min_coords.y=m->vnc_data.v[i].y;
    if (m->vnc_data.v[i].z < m->min_coords.z) m->min_coords.z=m->vnc_data.v[i].z;
  }

  _cleap_generate_edges_hash(m, dt);

  m->status = CLEAP_SUCCESS;
  m->wireframe = 0;
  m->solid = 1;

  return CLEAP_SUCCESS;
}

/*************************END CLEAP FUNCTIONS*********************/
int main(int argc, char* argv[])
{
  if (argc < 2) {
    fprintf(stderr, "usage: <exe>, <int: factor>\n");
    return -1;
  }

  int factor = atoi(argv[1]);
  if (factor < 0) {
    fprintf(stderr, "factor shouldn't be less than zero\n");
    exit(-1);
  }
  if (NUMPOINTS % factor != 0) {
    fprintf(stderr, "factor must be divisible by total number of points\n");
    exit(-1);
  }
  
  const int pop = NUMPOINTS / factor;
  Pair *DTarray = (Pair *)malloc(NUMPOINTS * sizeof(Pair));
  PointGenerator(DTarray);

  float *slabs = (float*)malloc((pop+3) * DIM * sizeof(float));

  qsort(DTarray, NUMPOINTS, sizeof(Pair), compareX);  
  slabs[2*pop] = -1; slabs[2*pop+1] = -1;
  slabs[2*pop+2] = 2; slabs[2*pop+3] = -1;
  slabs[2*pop+4] = .5; slabs[2*pop+5] = 2;

  DelaunayTriangulation dt;

  for(int i = 0; i < factor; i++) { 
    _cleap_mesh *m = new _cleap_mesh();
    if(cleap_init_no_render() != CLEAP_SUCCESS) {
      fprintf(stderr, "Failed to initialize correcly.\n");
      exit(-1);
    }

    slabPartition(DTarray, slabs, factor, i * pop, pop);    

    //fix this later
    dt.Initialize(-1, -1, 
		  2, -1, 
		  .5, 2);
    
    //create mesh
    for (int i = 0 ; i < pop; i++)
      dt.AddPoint(slabs[2*i], slabs[2*i+1]);

    //find faces of triangles in mesh
    dt.FindFaces(slabs, pop);
    
    //load cleap mesh
    if(load_mesh_host(pop, slabs, m, dt) != CLEAP_SUCCESS) {
      fprintf(stderr, "Failed to load mesh on host"); 
      exit(-1);
    }
    if(_cleap_device_load_mesh(m) != CLEAP_SUCCESS) {
      fprintf(stderr, "Failed to load mesh on device"); 
      exit(-1);
    } 
dt.VerifyResults(slabs, pop);
    //m = cleap_load_mesh("chair_0001.off");
    if(cleap_delaunay_transformation(m, CLEAP_MODE_2D) != CLEAP_SUCCESS) {
      fprintf(stderr, "Failed to run Delaunay"); 
      exit(-1);
    }

    if(i == 0) 
      cleap_save_mesh(m, "outMesh.off");
   
    dt.VerifyResults(slabs, pop); //error check
    //clear out vectors for next round
    dt.Clear();
    cleap_clear_mesh(m);
  } 

  //DT.WriteOutTriangle("kristi.vtk");

  free(slabs); 
  free(DTarray);

  return 0;
}
