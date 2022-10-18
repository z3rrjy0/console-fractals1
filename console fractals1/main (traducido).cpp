#include <Windows.h>
#include <iostream>
#include <vector>
#include <thread>
#include <map>
#include <sstream>
#include "CImg.h"
#include <complex>

constexpr signed char CUSTOM = -1;
constexpr signed char NORMAL = 1;
constexpr signed char MAP = 0;

constexpr long double pi = 3.141592653589793;
constexpr long double e = 2.718281828459045;

namespace ci = cimg_library;

unsigned int max_iterations = 255;

// COORD amplification = { 20.0f, 20.0f };

SHORT imageCount = 0;
bool scaleImage = true;
double long amplification;
COORD origin;
bool showCoordAndOrigin = false;
constexpr COORD topLeft = { 0, 0 };
CONSOLE_SCREEN_BUFFER_INFO csbi;
COORD fullscreen = { 480, 360 };
static const HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
std::map<std::string, SHORT> colours = {
  { "dblue", 1},
  { "dgreen", 2},
  { "lblue", 3},
  { "red", 4},
  { "lmagenta", 5},
  { "mango", 6},
  { "dwhite", 7},
  { "gray", 8},
  { "cian", 9},
  { "lgreen", 10},
  { "vintage", 11},
  { "salmon", 12},
  { "dmagenta", 13},
  { "meat", 14},
  { "white", 15}
};

void wait(const SHORT& s) {
  std::stringstream ss;
  ss << "timeout " << s << " > NULL";
  system(ss.str().c_str());
}

void getWindowSize() {
  if (!GetConsoleScreenBufferInfo(hOut, &csbi))
    // Gestió d'errors
    exit(-1); // abort();
}

COORD getHalfScreenSize() {
  COORD newCsbiDwSize;
  getWindowSize();
  newCsbiDwSize = csbi.dwSize;
  newCsbiDwSize.X /= 2;
  newCsbiDwSize.Y /= 2;
  return newCsbiDwSize;
}

template <typename T>
int size(T& v) {
  return sizeof(v) / sizeof(v[0]);
}

BOOL setCursorPosition(const int& x, const int& y) {
  std::cout.flush();
  const COORD c = { short(x), short(y) };
  return SetConsoleCursorPosition(hOut, c);
}

bool isVirtualKeyPressed(const int& virtKey) {
  return GetAsyncKeyState(virtKey) & 0x8000;
}

void changeColour(const int& colour) {
  SetConsoleTextAttribute(hOut, colour);
}

void drawSpot(const int& x, const int& y, const char& toDisplay, const int& colour = colours.at("white")) {
  if (x >= csbi.dwSize.X || x < 0 || y >= csbi.dwSize.Y || y < 0)
    return;
  setCursorPosition(x, y);
  changeColour(colour);
  std::cout << toDisplay;
  changeColour(colours.at("white"));
}

void cls() {
  std::cout.flush();
  getWindowSize();
  DWORD length = csbi.dwSize.X * csbi.dwSize.Y;
  DWORD written;
  FillConsoleOutputCharacter(hOut, TEXT(' '), length, topLeft, &written);
  FillConsoleOutputAttribute(hOut, csbi.wAttributes, length, topLeft, &written);
  SetConsoleCursorPosition(hOut, topLeft);
}

float mapNumber(const int& x, const int& in_min, const int& in_max, const int& out_min, const int& out_max) {
  try {
    return (x - in_min) * (out_max - out_min) / float(in_max - in_min + out_min);
  }
  catch (const std::exception& e) {
    cls();
    std::cout << "An error has ocurred!:\n\n" << e.what();
    exit(-1);
  }
}

float returnValueEquation(const std::vector<float>& equation, const float& x) {
  float toReturn = 0;
  unsigned int degree = 0;
  for (const float& monomyal : equation)
    toReturn += float(monomyal * pow(x, degree++));
  return toReturn;
}

float returnValueMapEquation(const std::map<float, float>& mapEquation, const float& x) {
  float toReturn = 0;
  for (const auto& monomyal : mapEquation)
    toReturn += monomyal.first * pow(x, monomyal.second);
  return toReturn; // *amplification); // Might be wrong
}

float returnValueCustomEquation(const float& x) { // Make slopeCustomEquation return the first derivative of this func
  return cos(x * 3);
  
  // return x % int(amplification); // (x ha de ser un valon sencer)
  float toReturn = 0;
  // Do something with toReturn
  return toReturn; // *amplification);
  return float(-2 * pow(x + 3, 2) * x - 5);
}

char charSlope(const float& slope) {
  // slope *= amplification;
  if (abs(slope) >= 17)  return '|';
  if (abs(slope) <= 3.5) return '-';
  if (slope > 0)         return '/';
                         return '\\';
}

char slopeEquation(std::vector<float> equation, const float& x) { // Still developing
  size_t eqS = equation.size();

  for (unsigned int i = 1; i < eqS;) // Derivar la funció
    equation[i - 1] = equation[i] * i++;
  equation.erase(equation.end() - 1);

  return charSlope(returnValueEquation(equation, x));
}

char slopeMapEquation(std::map<float, float> mapEquation, const float& x) { // Need an update

  for (auto it = mapEquation.begin(); it != mapEquation.end(); it++) {
    if (it->second != 0)
      mapEquation.erase(it);
    else
      it->second--;
  }
  /*
  unsigned int degree = 0; // May not even need the polynomial degree
  for (auto& element : mapEquation)
    // degree = (degree < element.second) ? degree : unsigned int(element.second);
    // if (degree < element.second) degree = element.second;
    degree = max(degree, element.second);

  // Unfinished
  */

  return charSlope(returnValueMapEquation(mapEquation, x));
}

char slopeCustomEquation(const float& x) { // Retorna el valor de la primera derivada (funció constant)
  return charSlope(float(-sin(x * 3) * amplification));
}

void plotGraph(const std::vector<float>& equation, const std::map<float, float>& mapEquation, const signed char& which) { // Can save the first derivative from the beginning
  LockWindowUpdate(GetConsoleWindow()); // Not sure if this actually locks window update
  // setCursorToPosition(0, csbi.dwSize.Y - 1);
  // std::cout << "\x1B[2J\x1B[H";
  cls(); // this func already calls getWindowSize();
  for (int i = 0; i < csbi.dwSize.Y; i++)
    drawSpot(origin.X, i, '|', colours.at("mango")); // The Y axis (|)
  float tmpResult;
  for (int i = 0; i < csbi.dwSize.X; i++) {
    drawSpot(i, origin.Y, '-', colours.at("mango")); // The X axis (-)
    char lineType = '*'; // Can remove the  = '*' part
    float tmpX = float((i - origin.X) / amplification);
    if (which == -1) { // CUSTOM
      tmpResult = float(returnValueCustomEquation(tmpX) * amplification);
      lineType = slopeCustomEquation(tmpX);
    }
    else if (which) { // NORMAL
      tmpResult = float(returnValueEquation(equation, tmpX) * amplification);
      lineType = slopeEquation(equation, tmpX);
    }
    else { // MAP
      tmpResult = float(returnValueMapEquation(mapEquation, tmpX) * amplification);
      lineType = slopeMapEquation(mapEquation, tmpX);
    }
    // tmpResult *= amplification / 2.5f;
    tmpResult /= 2.5f;

    drawSpot(i, origin.Y - int(tmpResult), lineType, colours.at("lblue")); // El gràfic (el guió pot ser un astirisc)
  }
  if (origin.X >= csbi.dwSize.X || origin.Y >= csbi.dwSize.Y) // Si l'origen es troba fora de la finestra, dibuixa l'x tan a prop com sigui possible
    drawSpot(max(min(origin.X, csbi.dwSize.X - 1), 1), max(min(origin.Y, csbi.dwSize.Y - 1), 1), 'X', colours.at("red"));
  else drawSpot(origin.X, origin.Y, 'O', colours.at("red")); // Sinó, dibuixa 'O'

  setCursorPosition(0, 0); // Print the amplification on the top left corner
  std::cout << "Amplification = " << amplification;
  LockWindowUpdate(NULL);
}

void promptData() {
  std::string toChange;
  std::stringstream newCoord;
  COORD half = getHalfScreenSize();
  cls();
  std::cout << "Write \"x\" for not changing\n";
  std::cout << "Previous amplification was " << amplification;
  std::cout << "\nAmplification = "; std::cin >> toChange;
  if (toChange != "x") amplification = stold(toChange);
  std::cout << "Previous origin in X was " << origin.X; // (half.X - origin.X) / amplification; // Amb la posició de la càmera // With the cam pos
  std::cout << "\nOrigin in X = "; std::cin >> toChange;
  if (toChange != "x") origin.X = stoi(toChange); // half.X - stoi(toChange) * amplification; // Amb la posició de la càmera // with the cam pos
  std::cout << "Previous origin in Y was " << origin.Y; // (half.Y - origin.Y) / amplification; // Amb la posició de la càmera // With the cam pos
  std::cout << "\nOrigin in Y = "; std::cin >> toChange;
  if (toChange != "x") origin.Y = stoi(toChange); // half.Y - stoi(toChange) * amplification; // Amb la posició de la càmera // With the cam pos
  std::cout << "Previous iteration limit was " << max_iterations;
  std::cout << "\nNew maximum iteration number: "; std::cin >> toChange;
  if (toChange != "x") max_iterations = stoi(toChange);
  std::cout << "Next fullscreen size in X (previous was " << fullscreen.X << "): "; std::cin >> toChange;
  if (toChange != "x") fullscreen.X = atoi(toChange.c_str());
  std::cout << "Next fullscreen size in Y (previous was " << fullscreen.Y << "): "; std::cin >> toChange;
  if (toChange != "x") fullscreen.Y = atoi(toChange.c_str());

  // origin = halfScreen - origin;
}

void hueToRGB(float H, float(&colour)[3]) {
  while (H > 360) H -= 360;
  while (H < 0) H += 360;

  float X = float(1 - abs(fmod(H / 60.0, 2) - 1)),
    r, g, b;

  if (H >= 0 && H < 60)         r = 1.0f, g = X, b = 0;
  else if (H >= 60 && H < 120)  r = X, g = 1.0f, b = 0;
  else if (H >= 120 && H < 180) r = 0, g = 1.0f, b = X;
  else if (H >= 180 && H < 240) r = 0, g = X, b = 1.0f;
  else if (H >= 240 && H < 300) r = X, g = 0, b = 1.0f;
  else                          r = 1.0f, g = 0, b = X;

  // const float colour2[] = { r * 255, g * 255, b * 255 };
  colour[0] = (r * 255) / 2;
  colour[1] = (g * 255) / 2;
  colour[2] = (b * 255);
  // return { r * 255, g * 255, b * 255 };
}

void HSVtoRGB(float HSV[3], float(&colour)[3]) { // Aquesta funció canvia la forma en la que l'ordinador interpreta un 'color'
  while (HSV[0] > 360) HSV[0] -= 360;
  while (HSV[0] < 0)   HSV[0] += 360;
  if (HSV[1]>100 || HSV[1] < 0 || HSV[2]>100 || HSV[2] < 0) return;

  float s = HSV[1] / 100;
  float v = HSV[2] / 100;
  float C = s * v;
  float X = float(C * (1 - abs(fmod(HSV[0]/ 60.0, 2) - 1)));
  float m = v - C;
  float r, g, b;

       if (HSV[0] >= 0   && HSV[0] < 60)  r = C, g = X, b = 0;
  else if (HSV[0] >= 60  && HSV[0] < 120) r = X, g = C, b = 0;
  else if (HSV[0] >= 120 && HSV[0] < 180) r = 0, g = C, b = X;
  else if (HSV[0] >= 180 && HSV[0] < 240) r = 0, g = X, b = C;
  else if (HSV[0] >= 240 && HSV[0] < 300) r = X, g = 0, b = C;
  else                                    r = C, g = 0, b = X;

  colour[0] = (r + m) * 255;
  colour[1] = (g + m) * 255;
  colour[2] = (b + m) * 255;
}

unsigned int mandelbrot(std::complex<long double>& c) {
  std::complex<long double> z = c;

  // Julia sets:

  // c = { -.79L, .15L };
  // c = { -.162L, 1.04L };
  // c = { .3, -.1 };
  // c = { -1.476, 0 };
  // c = { -.12, .77 };
  // c = { .28, .008 };
  // c = { -.7269, .1889 };

  for (unsigned int i = 0; i < max_iterations; ++i) {
    if (pow(z.real(), 2) + pow(z.imag(), 2) > 4.0) return i * 500; // mapNumber(i, 0, max_iterations, 1, 25); // i
    z = pow(z, 2) + c; // z = z^2 + c
  }
  return max_iterations*500;
}

std::complex<long double> newtonP(std::complex<long double>& z) {

  // return log(z);

  // return pow(z, e) * pow(e, z); // 10

  // return pow(z, 3) - z; // 9

  // return pow(pow(z, 4) - 1.L, 2) * pow(e, pow(z, 4)); // 8

  // return pow(z * e, pow(z, 6)); // 7

  // return (pow(z, 3) - 1.L) * pow(e, pow(z, 3)); // 6

  // return pow(z * e, pow(z, 3)); // 5 // slooow

  // return pow(z, 3) - z; // 4

  // return pow(z, 6) + pow(z, 3) - 1.L; // 3

  // return pow(z, 3) - 2.L * z + 2.L; // 2

  // return pow(z, 8) + 15.L * pow(z, 4) - 16.L; // 1

  // return pow(z, 15) - 1.L;

  return pow(z, 3) - 1.L;
  // z^3 - 1
}

std::complex<long double> newtonPPrime(std::complex<long double>& z) {

  // return 1.L / z;

  // return pow(e, z) * pow(z, e - 1.L) * (z + e); // 10

  // return 3.L * pow(z, 3) - 1.L; // 9

  // return 4.L * pow(e, pow(z, 4)) * pow(z, 3) * (pow(z, 8) - 1.L); // 8

  // return pow(e, pow(z, 6)) * pow(z, pow(z, 6) + 5.L) * (6.L * log10(z) + 7.L); // 7

  // return 3.L * pow(e, pow(z, 3)) * pow(z, 5); // 6

  // return pow(e, pow(z, 3)) * (3.L * pow(z, 3) + 1.L); // 5 // slooow

  // return 3.L * pow(z, 2) - 1.L; // 4

  // return 6.L * pow(z, 5), + 3.L * pow(z, 2); // 3

  // return 3.L * pow(z, 2) - 2.L; // 2

  // return 8.L * pow(z, 7) + 15.L * 4.L * pow(z, 3); // 1

  // return 15.L * pow(z, 14);

  return 3.L * pow(z, 2);
  // 3z^2 (p'(z)) = (z^3 - 1) dy/dx

  // z *= { 1/2., 0 }; // Generalització alternativa
}

unsigned int newton(const std::complex<long double>& c, float(&colour)[3]) { // Retorna un valor sencer
  std::complex<long double> z = c; // / amplification; // (c.X / amplification, c.Y / amplification);

  // std::complex<long double> p;

  float Tol = .0001f;
  std::complex<long double> r1(1, 0);
  std::complex<long double> r2(-0.5, sin(2 * pi / 3));
  std::complex<long double> r3(-0.5, -sin(2 * pi / 3));

  unsigned int i = 0;
  for (; i <= max_iterations && (abs(z - r1) >= Tol) && (abs(z - r2) >= Tol) && (abs(z - r3) >= Tol); ++i)
    z -= newtonP(z) / newtonPPrime(z); // z = z - P(z) / P'(z);

  colour[0] = 0;
  colour[1] = 0;
  colour[2] = 0;

  if (abs(z - r1) <  Tol) colour[0] = float(max_iterations - i);
  if (abs(z - r2) <= Tol) colour[1] = float(max_iterations - i);
  if (abs(z - r3) <= Tol) colour[2] = float(max_iterations - i);

  return i;
}

float collatz(const std::complex<long double>& c) {

  return float(((c*pow(cos(pi*c/2.L),2))/2.L+(3.L*c+1.L)+pow(sin(pi*c/2.L),2)/2.L).real());

  return float(((2.L + 7.L * c - (2.L + 5.L * c) * cos(pi * c)) / 4.L).real()); // Retorna la part real del nombre complex (2 + 7c - (2 + 5c)*cos(pi*c)) / 4

  // return float(((pow(-1, c) + 1.L) / 2.L * (c / 2.L) - (pow(-1, c) - 1.L) / 2.L * (3.L * c + 1.L)).real());
}

unsigned int burningShip(const std::complex<long double>& c) {
  std::complex<long double> z = c;

  for (unsigned int i = 0; i < max_iterations; i++) {
    
    if (abs(z) > 2) return i;

    // Zn + 1 = (| Zr | +i | Zi | ) ^ 2 + c
    z = pow(std::complex<long double>(abs(z.real()), abs(z.imag())), 2) + c;
  }

  return max_iterations;
}

unsigned int customF1(const std::complex<long double>& c) {
  std::complex<long double> z = c; // Initialize z as c

  for (unsigned int i = 0; i < max_iterations; ++i) {

    if (pow(abs(z), 2) > 4) return i; // mapNumber(i, 0, max_iterations, 1, 25); // i
    
    // z = (1.L-pow(z, 3)/6.L) / pow(z-pow(z, 2)/2.L, 2) + c;
    // z = (2.L + 7.L * z - (2.L + 5.L * z) * cos(pi * z)) / 4.L;
    // z -= a * ((pow(z, 4) - 1.L) / (4.L * pow(z, 3))); // Featherino
    // z = z - (pow(z, 2) - 1.L) / (a * (pow(z, 2) + 1.L)); // a was  2 (real)
    // z = conj(pow(z, 2)) + c; // Tricorn
    // z = conj(abs(pow(z, 2).real())) + c;
    // z = pow(sin(z), c) - pow(c, cos(z));
    z = pow(c, z) - pow(z, c);
  }
  return max_iterations;
}

unsigned int fractal(std::complex<long double>& c, float(&colour)[3]) { // Retorna el color del punt a on s'està calculant la funció
  c /= amplification;

  // ('colour' només s'aplica per a Zn+1 = z^3 + 1)
  // return newton(c, colour); //int toReturn = colour[0]; HSVtoRGB(colour, colour);
  // return toReturn;
  // (comentat però funciona)

  unsigned int tmpResult = mandelbrot(c); // També per als conjunts de Julia // Update
  // unsigned int tmpResult = collatz(c); // Doesn't plot what it's meant to do // Should change from "unsigned int" to "float"
  // unsigned int tmpResult = customF1(c);
  // unsigned int tmpResult = burningShip(c);

  hueToRGB(float(tmpResult * 360. / max_iterations), colour);
  return tmpResult;
}

const float colourR(const int& i) { // Blau
  return i / 3.f;
}

const float colourG(const int& i) { // Rosa
  return (255 + i) / 3.f;
}

const float colourB(const int& i) { // Groc
  return (510 + i) / 3.f;
}

void plotFractal(const bool saveAsImage) {
  if (saveAsImage) {
    goto afterBegin;
  begin:;
    // plotFractal(false); // cls(); si mostra percentatge
  afterBegin:;
    std::stringstream ss;
    ss << "mandelbrot" << rand() % 1000 + 9999 << rand() % 1000 + 9999 << rand() % 1000 + 9999 << ".bmp";
    std::cout << "Saving as " << ss.str();
    long double preAmp;
    COORD preOrigin;
    preAmp = amplification;
    preOrigin = origin;
    if (scaleImage) {
      amplification *= (fullscreen.X / csbi.dwSize.X + fullscreen.Y / csbi.dwSize.Y) / 2;
      origin.X *= fullscreen.X / csbi.dwSize.X;
      origin.Y *= fullscreen.Y / csbi.dwSize.Y;
    }
    // cls();
    float colour[3];
    ci::CImg<float> image(fullscreen.X, fullscreen.Y, 1, 3, 0);
    std::complex<long double> c;
    std::stringstream ss2;
    for (SHORT i = 0; i < fullscreen.Y; i++) {
      for (SHORT j = 0; j < fullscreen.X; j++) {
        c = { long double(j - origin.X), long double(i - origin.Y) };
        // max_iterations - fractal(C)
        fractal(c, colour); // mapNumber(isInCustomFractal(c, max_iterations), 0, 255, 0, 100);

#pragma region black-white-dblue
        /*
        float RG = (tmpResult * 2 < max_iterations) ? tmpResult * 2 : max_iterations - tmpResult*.94f*2.; // mapNumber(tmpResult, 0, max_iterations, 0, 240) * 2;
        float B = 0.;
        float colour[3];
        if (tmpResult * 2 < max_iterations) {
          colour[0] = RG;
          colour[1] = RG;
          colour[2] = RG;
        }
        else {
          if (tmpResult < max_iterations*.94f) { // 240
            colour[0] = max_iterations - tmpResult;
            colour[1] = colour[0];
            colour[2] = max_iterations;
          }
          else {
            tmpResult -= tmpResult*.94f; // mapNumber(tmpResult - max_iterations*.94f, 0, max_iterations*.06f, 0, 145); // 240 - 255 -> 0 - 15
            colour[0] = 0.;
            colour[1] = 0.;
            colour[2] = max_iterations - tmpResult;
          }
        }

        colour[0] = mapNumber(colour[0], 0, max_iterations, 0, 255);
        colour[1] = mapNumber(colour[1], 0, max_iterations, 0, 255);
        colour[2] = mapNumber(colour[2], 0, max_iterations, 0, 255);
        */
#pragma endregion

        // hueToRGB(tmpResult, colour);
        // const float colour[] = { RG, RG, 360. }; // { hueToRGB(tmpResult, 0), hueToRGB(tmpResult, 1), hueToRGB(tmpResult, -1) }; // { colourR(tmpResult), colourG(tmpResult), colourB(tmpResult) };
        // const float colour[3] = { tmpResult, tmpResult, tmpResult }; // temp
        // hueToRGB(tmpResult * 360. / max_iterations, colour);
        image.draw_point(j, i, colour); // Colour B, Pink, Y // Actually RGB
      }
      ss2.str(""); ss2 << i * 100.0 / fullscreen.Y << "%"; setCursorPosition(0, 4); std::cout << "          "; setCursorPosition(0, 4); std::cout << ss2.str();
    }
    float red[3] = { 255., 0., 0. };
    // Amp and if were here
    if (showCoordAndOrigin) {
      image.draw_line(origin.X, 0, origin.X, fullscreen.Y, red); // eix Y
      image.draw_line(0, origin.Y, fullscreen.X, origin.Y, red); // eix X
      if (origin.X >= fullscreen.X || origin.Y >= fullscreen.Y) image.draw_circle(max(min(origin.X, fullscreen.X - 1), 1), max(min(origin.Y, fullscreen.Y - 1), 1), 3, red);
      else image.draw_circle(origin.X, origin.Y, 3, red); // Dues linies pel punt central
    }
    amplification = preAmp;
    if (scaleImage)
      origin = preOrigin;
    image.save(ss.str().c_str());
    ci::CImgDisplay local(image, ss.str().c_str(), 0); // add a true argument for fullscreen mode
    local.move(fullscreen.X/2, fullscreen.Y/2);
    // while(!isVirtualKeyPressed(VK_RETURN)) local.wait(); return;
    int sens = 50;

    const long double sensAmp = 2.0f;
    long double preAmp2 = amplification;
    COORD halfScreen;
  before:; // Configuració de la finestra flotant
         if (isVirtualKeyPressed(VK_LEFT))  origin.X += sens; // UP
    else if (isVirtualKeyPressed(VK_RIGHT)) origin.X -= sens;
    else if (isVirtualKeyPressed(VK_UP))    origin.Y += SHORT(sens * 2.2f); // Sense 2.2f
    else if (isVirtualKeyPressed(VK_DOWN))  origin.Y -= SHORT(sens * 2.2f); // Sense 2.2f
    else if (isVirtualKeyPressed(VK_OEM_PLUS))  amplification *= sensAmp; // amplification = max(min(amplification * sensAmp, INT_MAX), 1);
    else if (isVirtualKeyPressed(VK_OEM_MINUS)) amplification /= sensAmp; // amplification = max(min(amplification / sensAmp, INT_MAX), 1);
    else if (isVirtualKeyPressed(VK_F5)) { if (scaleImage) { getWindowSize(); origin = csbi.dwSize; } else origin = fullscreen; origin.X /= 2; origin.Y /= 2; amplification = 20; }
    else if (isVirtualKeyPressed(VK_F6)) promptData();
    else if (isVirtualKeyPressed(VK_F11)) max_iterations--;
    else if (isVirtualKeyPressed(VK_F12)) max_iterations++;
    else if (isVirtualKeyPressed(VK_RETURN)) return;
    else goto before;
    if (preAmp2 != amplification) {
      halfScreen = getHalfScreenSize();
      origin.X = SHORT(halfScreen.X - (halfScreen.X - origin.X) * amplification / preAmp2);
      origin.Y = SHORT(halfScreen.Y - (halfScreen.Y - origin.Y) * amplification / preAmp2);
      preAmp2 = amplification;
    }
    goto begin;
    // plotFractal(fractal, mapFractal, which, max_iterations, saveAsImage);
    return;
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  LockWindowUpdate(GetConsoleWindow());
  // setCursorToPosition(0, csbi.dwSize.Y - 1);
  // std::cout << "\x1B[2J\x1B[H";
  cls();
  std::complex<long double> c;
  unsigned int tmpResult;
  float colour[3];

  ci::CImg<float> image(csbi.dwSize.X, csbi.dwSize.Y, 1, 3, 0);

  setCursorPosition(0, 0);
  for (SHORT i = 0; i < csbi.dwSize.Y; i++) {
    for (SHORT j = 0; j < csbi.dwSize.X; j++) {

      c = { long double(j - origin.X), (i - origin.Y) * 2.2f };
      tmpResult = fractal(c, colour);
      
      // tmpResult = mapNumber(tmpResult, 0, 1, 0, 255);
      /*
      stringstream ss;
      ss << "tmpResult = " << tmpResult;
      SetCursorPos(0, 10);
      changeColour(colours.at("white"));
      std::cout << ss.str();
      */ // Not the best part of the code to std::cout anything...

      changeColour(tmpResult);
      std::cout << char(219); // "█"
      // drawSpot(i, j, '*', tmpResult); // The fractal (*) // tmpResult -> colours.at("dblue")

      if (saveAsImage) { // Quan es guarda en una imatge, no es guarden els eixos de coordenades
        // colour = { colourR(tmpResult), colourG(tmpResult), colourB(tmpResult) };
        image.draw_point(j, i, colour);
      }

    }
    if (i != csbi.dwSize.Y - 1) std::cout << "\n";
  }
  if (saveAsImage) { // Hauria de moure's després de mstrar l'original
    std::stringstream imageFileName;
    image.save("file.bmp"); goto after;
    imageFileName << "file" << ++imageCount << ".bmp";
    image.save(imageFileName.str().c_str());
  after:;
    ci::CImgDisplay local(image, "Mandelbrot", 0, true);
    local.wait();
    std::cout << "Finished!!!\n";
  }
  if (showCoordAndOrigin) {
    for (int i = 0; i < csbi.dwSize.Y; i++) drawSpot(origin.X, i, '|', colours.at("mango")); // The Y axis (|)
    for (int i = 0; i < csbi.dwSize.X; i++) drawSpot(i, origin.Y, '-', colours.at("mango")); // The X axis (-)
    if (origin.X >= csbi.dwSize.X || origin.Y >= csbi.dwSize.Y) drawSpot(max(min(origin.X, csbi.dwSize.X - 1), 1), max(min(origin.Y, csbi.dwSize.Y - 1), 1), 'X', colours.at("red"));
    else drawSpot(origin.X, origin.Y, 'O', colours.at("red"));
  }
  setCursorPosition(0, 0);
  changeColour(colours.at("white"));
  std::cout << "Amplification = " << amplification;
  std::cout << "\nOrigin = " << origin.X << ", " << origin.Y;
  COORD half = getHalfScreenSize();
  half.X -= origin.X; half.Y -= origin.Y; // No operator "-=" between "COORD" and "COORD"
  std::cout << "\nCamera position = " << half.X / amplification << ", " << -half.Y / amplification << "\n";
  LockWindowUpdate(NULL);
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}

void printEq(const std::vector<float>& equation) { // temp
  for (size_t i = equation.size() - 1; i >= 0;i--) {
    if (equation[i] != 0) {
      std::cout << "\t";
      if (equation[i] != 1)
        std::cout << equation[i];
      if(i!=0)
        std::cout << "x" << i;
    }
  }
  std::cout << "\t<--\n";
}

int main() {
  
  // getWindowSize();
  // std::cout << csbi.dwSize.X << "\t" << csbi.dwSize.Y << "\n";
  // return 0;

  srand(unsigned int(time(NULL)));

  bool useFr = true; // true  -> Fractal
                      // false -> Gràfic

  bool saveAsImage = false; // Encara només per a fractals

  const signed char whichEq = CUSTOM; // CUSTOM / NORMAL / MAP

  std::vector<float> equation = { 0, 0, 1 };
  std::map<float, float> mapEquation = {
    {1, 3}
  };

  goto skipTest;
  /////////////////////////// startdelete
  printEq(equation);
  for (unsigned int i = 1; i < equation.size();)
    equation[i-1] = equation[i] * i++;
  equation.erase(equation.end() - 1);
  // equation[degree] = 0; // instead of the last line of code
  printEq(equation);
  std::cout << "\nFinished\n";
  return 0;
  /////////////////////////// enddelete
skipTest:;

  amplification = 20;
  
  cls(); // getWindowSize();
  origin = csbi.dwSize;
  origin.X /= 2; // Posició inicial de l'origen de coordenades
  origin.Y /= 2; // El centre de la consola de l'ordinador
  
  // origin.X = fullscreen.X / 2; // La posició inicial de l'origen de
  // origin.Y = fullscreen.Y / 2; // coordenades es el centre de la consola

  int sens; // Indica la sensibilitat de la fletxa
  long double sensAmp; // Indica la sensibilitat dels controls d'amplificació

  if (useFr) {
    plotFractal(saveAsImage); // , int(amplification));
    sens = 50;
    sensAmp = 1.5;
  }
  else {
    plotGraph(equation, mapEquation, whichEq);
    sens = csbi.dwSize.X * csbi.dwSize.Y * 3 / 3782; // 3
    sensAmp = 1.1;
  }

  long double preAmp = amplification;
  COORD halfScreen;

  while (!isVirtualKeyPressed(VK_ESCAPE)) { // Configuració de la consola
    if (isVirtualKeyPressed(VK_LEFT)) origin.X += sens; // UP (no la tecla)
    else if (isVirtualKeyPressed(VK_RIGHT)) origin.X -= sens;
    else if (isVirtualKeyPressed(VK_UP))    origin.Y += SHORT(sens / 2.2f); // Els caràcters de la consola no son quadrats, sinó rectangles,
    else if (isVirtualKeyPressed(VK_DOWN))  origin.Y -= SHORT(sens / 2.2f); // això evita que el dibuix surti estirat cam amunt
    else if (isVirtualKeyPressed(VK_OEM_PLUS))  amplification *= sensAmp; // amplification = max(min(amplification * sensAmp, INT_MAX), 1); // Min not one
    else if (isVirtualKeyPressed(VK_OEM_MINUS)) amplification /= sensAmp; // amplification = max(min(amplification / sensAmp, INT_MAX), 1); // Min not one
    else if (isVirtualKeyPressed(VK_F5)) { getWindowSize(); origin = csbi.dwSize; origin.X /= 2; origin.Y /= 2; amplification = 20; } // Restablir posició i amplificació
    else if (isVirtualKeyPressed(VK_F6)) promptData();
    else if (isVirtualKeyPressed(VK_F1) && useFr) saveAsImage = true;
    else if (isVirtualKeyPressed(VK_F11)) max_iterations--;
    else if (isVirtualKeyPressed(VK_F12)) max_iterations++;
    else continue; // Si es compleix algun condicional continua el codi, sinó torna al bucle
    // getWindowSize();
    if (preAmp != amplification) {
      halfScreen = getHalfScreenSize();
      origin.X = SHORT(halfScreen.X - (halfScreen.X - origin.X) * amplification / preAmp);
      origin.Y = SHORT(halfScreen.Y - (halfScreen.Y - origin.Y) * amplification / preAmp);
      preAmp = amplification;
    }
    if (useFr) {
      plotFractal(saveAsImage);
      saveAsImage = false;
    }
    else {
      plotGraph(equation, mapEquation, whichEq);
      sens = csbi.dwSize.X * csbi.dwSize.Y * 3 / 3782; // 3
    }
  }

  setCursorPosition(0, csbi.dwSize.Y - 1);
  return 0;
}
