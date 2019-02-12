////////////////////////////////////////////////////////////////////////
//
//   Source code based on asst.cpp by
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////

#pragma
#include "pch.h"
#include <vector>
#include <string>
#include <memory>
#include <stdexcept>

#include <GL/glew.h>
#include <GL/glut.h>

#include "cvec.h"
#include "matrix4.h"
#include "rigtform.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"

using namespace std;      // for string, vector, iostream, and other standard C++ stuff

// G L O B A L S ///////////////////////////////////////////////////

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = -2.0;      // y coordinate of the ground
static const float g_groundSize = 10.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;


struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;
  GLint h_uUseTexture;
  GLint h_uTexUnit;
  GLint h_aTexCoord;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;
  GLint h_aTexCoord0;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    //h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");
	h_uUseTexture = safe_glGetUniformLocation(h,"uUseTexture");
    h_uTexUnit = safe_glGetUniformLocation(h, "uTexUnit");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");
    h_aTexCoord = safe_glGetAttribLocation(h, "aTexCoord");

    checkGlErrors();
  }

};

static const int g_numShaders = 2;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic.vshader", "./shaders/diffuse.fshader"},
  {"./shaders/basic.vshader", "./shaders/specular.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

static shared_ptr<GlTexture> g_tex; // our global texture instance

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position, normal, and texture
struct VertexPNT {
  Cvec3f p, n;
  Cvec2f t;

  VertexPNT() {}
  VertexPNT(float x, float y, float z,
           float nx, float ny, float nz,
		   float tx, float ty)
    : p(x,y,z), n(nx, ny, nz), t(tx,ty)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPNT(const GenericVertex& v) {
    *this = v;
  }

  VertexPNT& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
	t = v.tex;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPNT *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPNT) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);
    safe_glEnableVertexAttribArray(curSS.h_aTexCoord);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPNT), FIELD_OFFSET(VertexPNT, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPNT), FIELD_OFFSET(VertexPNT, n));
	safe_glVertexAttribPointer(curSS.h_aTexCoord, 2, GL_FLOAT, GL_FALSE, sizeof(VertexPNT), FIELD_OFFSET(VertexPNT, t));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
    safe_glDisableVertexAttribArray(curSS.h_aTexCoord);
  }
};


// Vertex buffer and index buffer associated with the ground and cube geometry amd sphere geometry
static shared_ptr<Geometry> g_ground, g_cube, g_sphere;

// --------- Scene

static const Cvec3 g_light(0.0, 0.25, 4.0);    // define a light positions in world space//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
static Matrix4 g_skyRbt = Matrix4::makeTranslation(Cvec3(0.0, 0.25, 4.0));
static Matrix4 g_objectRbt[3] = {Matrix4::makeTranslation(Cvec3(-1,0,0)),Matrix4::makeTranslation(Cvec3(1,0,0)), Matrix4::makeTranslation(g_light)};  // Two cubes and a sphere
static Cvec3f g_objectColors[3] = {Cvec3f(1, 0, 0),Cvec3f(0, 0, 1),Cvec3f(1, 1, 0) };
static int g_activeCube = 0;

///////////////// END OF G L O B A L S //////////////////////////////////////////////////


static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPNT vtx[4] = {
    VertexPNT(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0,0,0),
    VertexPNT(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0,0,1),
    VertexPNT( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0,1,1),
    VertexPNT( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0,1,0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static void initCubes() {
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPNT> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  g_cube.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

static void initSphere() {
	int ibLen, vbLen;
	getSphereVbIbLen(8,8,vbLen, ibLen);

	// Temporary storage for sphere geometry
	vector<VertexPNT> vtx(vbLen);
	vector<unsigned short> idx(ibLen);

	makeSphere(0.5, 8, 8, vtx.begin(), idx.begin());
	g_sphere.reset(new Geometry(&vtx[0], &idx[0], vbLen, ibLen));
}

// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS175_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}

static void drawStuff() {
  // short hand for current shader state
  const ShaderState& curSS = *g_shaderStates[g_activeShader];

  // build & send proj. matrix to vshader
  const Matrix4 projmat = makeProjectionMatrix();
  sendProjectionMatrix(curSS, projmat);

  // use the skyRbt as the eyeRbt
  const Matrix4 eyeRbt =g_skyRbt;
  const Matrix4 invEyeRbt = inv(eyeRbt);

  const Cvec3 eyeLight = Cvec3(invEyeRbt * Cvec4(g_light, 1)); // g_light position in eye coordinates
  safe_glUniform3f(curSS.h_uLight, eyeLight[0], eyeLight[1], eyeLight[2]);
  safe_glUniform1i(curSS.h_uUseTexture,0); // Turn off textures for the ground and the first cube
  
  // draw ground
  // ===========
  //
  const Matrix4 groundRbt = Matrix4();  // identity
  Matrix4 MVM = invEyeRbt * groundRbt;
  Matrix4 NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
  g_ground->draw(curSS);

  // draw cubes
  // ==========
  MVM = invEyeRbt * g_objectRbt[0];
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_objectColors[0][0], g_objectColors[0][1], g_objectColors[0][2]);
  g_cube->draw(curSS);

  MVM = invEyeRbt * g_objectRbt[1];
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform1i(curSS.h_uUseTexture,1); // Turn on textures for the second cube
  g_cube->draw(curSS);

  // draw a sphere
  // ==========
  safe_glUniform1i(curSS.h_uUseTexture, 0); // Turn off textures for the sphere
  MVM = invEyeRbt * g_objectRbt[2];
  NMVM = normalMatrix(MVM);
  sendModelViewNormalMatrix(curSS, MVM, NMVM);
  safe_glUniform3f(curSS.h_uColor, g_objectColors[2][0], g_objectColors[2][1], g_objectColors[2][2]);
  g_sphere->draw(curSS);
}

static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}
static void motion(const int x, const int y) {
	const double dx = x - g_mouseClickX;
	const double dy = g_windowHeight - y - 1 - g_mouseClickY;

	if (g_mouseClickDown) {
		Matrix4 A = transFact(g_objectRbt[g_activeCube]) * linFact(g_skyRbt);
		Matrix4 M;

		if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
			// Do M to O with respect to A
			M = Matrix4::makeXRotation(-dy) * Matrix4::makeYRotation(dx);
			g_objectRbt[g_activeCube] = A * M*inv(A)*g_objectRbt[g_activeCube];
		}
		else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
			M = Matrix4::makeTranslation(Cvec3(dx, dy, 0) * 0.01);
			g_objectRbt[g_activeCube] = A * M*inv(A)*g_objectRbt[g_activeCube];
		}
		else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
			M = Matrix4::makeTranslation(Cvec3(0, 0, -dy) * 0.01);
			g_objectRbt[g_activeCube] = A * M*inv(A)*g_objectRbt[g_activeCube];
		}
		if (g_activeCube == 2) {
			g_light = g_objectRbt[g_activeCube]
		}
		else {}
		glutPostRedisplay(); // we always redraw if we changed the scene

	}

	g_mouseClickX = x;
	g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;
}

static void keyboard(const unsigned char key, const int x, const int y) {
  switch (key) {
  case 27:
    exit(0);                                  // ESC
  case 'h':
    cout << " ============== H E L P ==============\n\n"
    << "h\t\thelp menu\n"
    << "s\t\tsave screenshot\n"
    << "f\t\tToggle flat shading on/off.\n"
    << "o\t\tCycle object to edit\n"
    << "1\t\tDiffuse only\n"
	<< "2\t\tDiffuse and specular\n"
    << "drag left mouse to rotate\n" << endl;
    break;
  case 's':
    glFlush();
    writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
    break;
  case 'o':
	  if (g_activeCube == 0) g_activeCube = 1;
	  else if (g_activeCube == 1) g_activeCube = 2;
	  else g_activeCube = 0;
	  break;
  case '1':
	  g_activeShader = 0;
	  break;
  case '2':
	  g_activeShader = 1;
	  break;
  }
  glutPostRedisplay();
}

static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Lighting");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() {
  initGround();
  initCubes();
  initSphere();
}

static void loadTexture(GLuint texHandle, const char *ppmFilename) {
  int texWidth, texHeight;
  vector<PackedPixel> pixData;

  ppmRead(ppmFilename, texWidth, texHeight, pixData);

  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texHandle);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_SRGB, texWidth, texHeight,
               0, GL_RGB, GL_UNSIGNED_BYTE, &pixData[0]);
  checkGlErrors();
}


static void initTextures() {
  g_tex.reset(new GlTexture());
  
  loadTexture(*g_tex, "smiley.ppm");
 
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, *g_tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
}


int main(int argc, char * argv[]) {
  try {
    initGlutState(argc,argv);

    glewInit(); // load the OpenGL extensions

    if (!GLEW_VERSION_3_0)
      throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
    
    initGLState();
    initShaders();
    initGeometry();
    initTextures();

    glutMainLoop();
    return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
