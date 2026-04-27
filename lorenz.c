#include <math.h>
#include <stdlib.h>
#include <SDL2/SDL.h>
#include <time.h>

#define WIDTH 900
#define HEIGHT 600

#define COLOR_WHITE 0xffffff
#define h 0.01
#define STEPS 120000
#define NUMBER_OF_EQS 3
#define PROYECTION_PARAM 35
#define SCALE 15

#define OFFSET_X (WIDTH/2)
#define OFFSET_Y ((1.5)* HEIGHT/2)

int POINT_SIZE = 1;


typedef struct{
    double x, y, z;
}Point;

// Chen-Lee attractor
/*void f(double t, double y[], double dydt[]){

    double a = 5.0;
    double b = -10.0;
    double c = 1.0/3.0;
    double d = - 0.38;  

    dydt[0] = a * y[0] - y[1] * y[2];
    dydt[1] = b * y[1] + y[0] * y[2];
    dydt[2] = d * y[2] + y[1] * y[2] * c;
}
*/

// Lorenz atracttor
/*
void f(double t, double y[], double dydt[]){
    double sigma = 10.0;
    double rho   = 28.0;
    double beta  = 8.0/3.0;

    dydt[0] = sigma * (y[1] - y[0]);                 // dx/dt
    dydt[1] = y[0] * (rho - y[2]) - y[1];            // dy/dt
    dydt[2] = y[0] * y[1] - beta * y[2];             // dz/dt
}
*/

// Dadras Attractor
/*
void f(double t, double y[], double dydt[]){

    double a = 3.0;
    double b = 2.7;
    double c = 1.7;
    double d = 2.0;
    double e = 9.0; 

    double x = y[0];
    double yv = y[1];
    double z = y[2];

    dydt[0] = yv - a * x + b * yv * z;
    dydt[1] = c * yv - x * z + z;
    dydt[2] = d * x * yv - e * z;
}
*/
// Rösler attractor

void f(double t, double y[], double dydt[]){
    double a = 0.2;
    double b = 0.2;
    double c = 5.7;

    dydt[0] = -y[1] - y[2];            // dx/dt
    dydt[1] = y[0] + a * y[1];         // dy/dt
    dydt[2] = b + y[2] * (y[0] - c);   // dz/dt
}


Point* Runge_Kutta(double t, double y[]){
    Point* path = malloc(STEPS *sizeof(Point));
    double k1[NUMBER_OF_EQS], k2[NUMBER_OF_EQS], k3[NUMBER_OF_EQS], k4[NUMBER_OF_EQS], yt[NUMBER_OF_EQS];

    for(int step = 0; step < STEPS; step++){

        f(t, y, k1);

        for (int i = 0; i < NUMBER_OF_EQS; i++)
            yt[i] = y[i] + 0.5 * h * k1[i];
        f(t + 0.5 * h, yt, k2);

        for (int i = 0; i < NUMBER_OF_EQS; i++)
            yt[i] = y[i] + 0.5 * h * k2[i];
        f(t + 0.5 * h, yt, k3);

        for (int i = 0; i < NUMBER_OF_EQS; i++)
            yt[i] = y[i] + h * k3[i];
        f(t + h, yt, k4);

        for (int i = 0; i < NUMBER_OF_EQS; i++)
            y[i] += (h / 6.0) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
        t += h;
        path[step].x = y[0];
        path[step].y = y[1];
        path[step].z = y[2];
    }
    return path;
}

void apply_rotation(Point* point, double phi){
    double temp_point[3] = {point->x, point->y, point->z};
    double result[3] = {0, 0, 0};
    double rotation_matrix[3][3] = {
        {1, 0, 0},
        {0, cos(phi), -sin(phi)},
        {0, sin(phi),  cos(phi)}
    };
    double sum = 0;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            sum += rotation_matrix[i][j] * temp_point[j];
        }
        result[i] = sum;
        sum = 0;
    }

    point->x = result[0];
    point->y = result[1];
    point->z = result[2];
}


void draw_point(SDL_Surface *psurface, int x, int y){
    SDL_Rect rect = (SDL_Rect){x, y, POINT_SIZE, POINT_SIZE};
    SDL_FillRect(psurface, &rect, COLOR_WHITE);
}

void draw_point_3d(SDL_Surface *psurface, Point point){
    double dist = PROYECTION_PARAM;
    double factor = dist / (dist + point.z);
    int x_2d = point.x * factor * SCALE + OFFSET_X;
    int y_2d = point.y * factor * SCALE + OFFSET_Y;

    draw_point(psurface, x_2d, y_2d);
}


Point* Rotate_Solution(Point* path, double angle){
    for(int i = 0; i < STEPS; i++){
        apply_rotation(&path[i], angle);
    }
    return path;
}

void print_path(SDL_Surface *psurface, Point* path, double phi, int limit){
    if(limit < STEPS){
        for(int i = 0; i < limit; i++){
            Point p = path[i];
            apply_rotation(&p, phi);
            draw_point_3d(psurface, p);
        }
    }else{
        for(int i = 0; i < STEPS; i++){
            Point p = path[i];
            apply_rotation(&p, phi);
            draw_point_3d(psurface, p);
        }
    }
}

int main(){
    double t = 0.0;
    double phi = 2.0;
    /* Dadras/Rössler/Lorenz Initial Conditions
    double y[NUMBER_OF_EQS] = {1.0, 1.0, 1.0};
    */
    double y[NUMBER_OF_EQS] = {1.0, 1.0, 1.0};
    Point* path = Runge_Kutta(t, y);
    
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        printf("SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow("Attractors",
                                            SDL_WINDOWPOS_CENTERED, 
                                            SDL_WINDOWPOS_CENTERED, 
                                            WIDTH, 
                                            HEIGHT, 
                                            0);

    SDL_Surface *psurface = SDL_GetWindowSurface(window);
    //init_colors(psurface);
    if (!window) {
        printf("SDL_CreateWindow error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    int running = 1;
    SDL_Event e;
    int limit = 1;
    while (running) {
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT)
                running = 0;
        }
        SDL_FillRect(psurface, NULL, 0x000000);
        print_path(psurface, path, phi, limit);
        SDL_UpdateWindowSurface(window);
        // phi += 0.1;
        limit++;
        
    }

    SDL_DestroyWindow(window);
    SDL_Quit();
    free(path);

    return 0;
}