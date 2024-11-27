/// 
//  mandel.c
//  Based on example code found here:
//  https://users.cs.fiu.edu/~cpoellab/teaching/cop4610_fall22/project3.html
//
//  Converted to use jpg instead of BMP and other minor changes
//  
///

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <semaphore.h>
#include <fcntl.h>
#include <math.h>
#include "jpegrw.h"

// Function prototypes
int iterations_at_point(double x, double y, int max);   // Calculate iterations for a given point
int iteration_to_color(int i, int max);                // Map iteration count to color
void compute_image(imgRawImage *img, double xmin, double xmax, double ymin, double ymax, int max); // Compute Mandelbrot fractal
void show_help();                                      // Display help message
void generate_movie(double xcenter, double ycenter, double xscale, int image_width, 
                    int image_height, int max, int num_frames, int num_children); // Generate a Mandelbrot zooming movie

int main(int argc, char *argv[]) {
    const char *outfile = "mandel.jpg"; // Default output file for single image
    double xcenter = -0.5, ycenter = 0; // Default center coordinates of the fractal
    double xscale = 4, yscale = 0;     // Default scale and aspect ratio
    int image_width = 1000, image_height = 1000, max = 1000; // Default resolution and max iterations

    // Check if "--movie" mode is enabled
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--movie") == 0) {
            int num_frames = 50, num_children = 4; // Default number of frames and processes
            char c;

            // Parse movie-specific options
            while ((c = getopt(argc, argv, "x:y:s:W:H:m:f:p:h")) != -1) {
                switch (c) {
                    case 'x': xcenter = atof(optarg); break;      // Set X center
                    case 'y': ycenter = atof(optarg); break;      // Set Y center
                    case 's': xscale = atof(optarg); break;       // Set scale
                    case 'W': image_width = atoi(optarg); break;  // Set image width
                    case 'H': image_height = atoi(optarg); break; // Set image height
                    case 'm': max = atoi(optarg); break;          // Set max iterations
                    case 'f': num_frames = atoi(optarg); break;   // Set number of frames
                    case 'p': num_children = atoi(optarg); break; // Set number of processes
                    case 'h': show_help(); return 0;              // Display help and exit
                }
            }
            generate_movie(xcenter, ycenter, xscale, image_width, image_height, max, num_frames, num_children);
            return 0; // Exit after generating the movie
        }
    }

    // Parse single-image options
    char c;
    while ((c = getopt(argc, argv, "x:y:s:W:H:m:o:h")) != -1) {
        switch (c) {
            case 'x': xcenter = atof(optarg); break;      // Set X center
            case 'y': ycenter = atof(optarg); break;      // Set Y center
            case 's': xscale = atof(optarg); break;       // Set scale
            case 'W': image_width = atoi(optarg); break;  // Set image width
            case 'H': image_height = atoi(optarg); break; // Set image height
            case 'm': max = atoi(optarg); break;          // Set max iterations
            case 'o': outfile = optarg; break;            // Set output file
            case 'h': show_help(); return 0;              // Display help and exit
        }
    }

    // Compute the aspect ratio
    yscale = xscale / image_width * image_height;

    // Print the rendering parameters
    printf("Rendering Mandelbrot image with parameters:\n");
    printf("Center: (%lf, %lf), Scale: %lf, Resolution: %dx%d, Max Iterations: %d\n", 
           xcenter, ycenter, xscale, image_width, image_height, max);

    // Initialize the image
    imgRawImage *img = initRawImage(image_width, image_height);

    // Compute the Mandelbrot fractal
    compute_image(img, xcenter - xscale / 2, xcenter + xscale / 2, 
                  ycenter - yscale / 2, ycenter + yscale / 2, max);

    // Save the fractal to a file
    storeJpegImageFile(img, outfile);

    // Free allocated memory
    freeRawImage(img);
    printf("Image saved to %s\n", outfile);

    return 0;
}

// Generate a zooming Mandelbrot movie
void generate_movie(double xcenter, double ycenter, double xscale, int image_width, 
                    int image_height, int max, int num_frames, int num_children) {
    // Remove any existing semaphore with the same name
    sem_unlink("/mandel_sem");

    // Create a new semaphore
    sem_t *sem = sem_open("/mandel_sem", O_CREAT | O_EXCL, 0644, num_children);
    if (sem == SEM_FAILED) {
        perror("sem_open failed");
        exit(1);
    }

    double zoom_factor = 0.90; // Zoom factor for each frame

    // Generate each frame
    for (int frame = 0; frame < num_frames; frame++) {
        sem_wait(sem); // Wait for a free slot

        // Fork a new child process
        if (fork() == 0) {
            char filename[256];
            snprintf(filename, sizeof(filename), "frame%02d.jpg", frame);

            // Calculate the bounds for the current frame
            double scale = xscale * pow(zoom_factor, frame);
            double xmin = xcenter - scale / 2;
            double xmax = xcenter + scale / 2;
            double ymin = ycenter - scale / 2;
            double ymax = ycenter + scale / 2;

            // Debugging: Log frame details
            printf("Frame %d: xmin=%f, xmax=%f, ymin=%f, ymax=%f, scale=%f\n", 
                   frame, xmin, xmax, ymin, ymax, scale);

            // Generate the fractal image for the frame
            imgRawImage *img = initRawImage(image_width, image_height);
            compute_image(img, xmin, xmax, ymin, ymax, max);
            if (storeJpegImageFile(img, filename) != 0) {
                fprintf(stderr, "Failed to save frame: %s\n", filename);
            }
            freeRawImage(img);

            // Release the semaphore slot and exit child process
            sem_post(sem);
            exit(0);
        }
    }

    // Wait for all child processes to finish
    while (wait(NULL) > 0);

    // Cleanup the semaphore
    sem_close(sem);
    sem_unlink("/mandel_sem");

    // Create a movie from the frames
    printf("All frames generated! Creating movie...\n");
    system("ffmpeg -y -framerate 24 -i frame%02d.jpg mandelmovie.mp4");
    printf("Movie created: mandelmovie.mp4\n");
}

// Compute the Mandelbrot fractal
void compute_image(imgRawImage *img, double xmin, double xmax, double ymin, double ymax, int max) {
    int width = img->width;
    int height = img->height;

    // Iterate through each pixel in the image
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double cx = xmin + x * (xmax - xmin) / (width - 1);
            double cy = ymin + y * (ymax - ymin) / (height - 1);
            int iterations = iterations_at_point(cx, cy, max); // Calculate iterations
            int color = iteration_to_color(iterations, max);  // Map iterations to color
            setPixelCOLOR(img, x, y, color);                 // Set pixel color
        }
    }
}

// Calculate the number of iterations at a given point
int iterations_at_point(double x, double y, int max) {
    int count = 0;
    double zx = 0.0, zy = 0.0;
    while (zx * zx + zy * zy < 4.0 && count < max) {
        double temp = zx * zx - zy * zy + x;
        zy = 2.0 * zx * zy + y;
        zx = temp;
        count++;
    }
    return count;
}

// Map the iteration count to a color
int iteration_to_color(int i, int max) {
    if (i == max) return 0x000000; // Black for points inside the Mandelbrot set
    int r = (i * 9) % 256; // Red gradient
    int g = (i * 15) % 256; // Green gradient
    int b = (i * 25) % 256; // Blue gradient
    return (r << 16) | (g << 8) | b; // Combine RGB
}

// Display help message
void show_help() {
    printf("Usage: mandel [options]\n");
    printf("Options:\n");
    printf("  -x <double>   X center (default=-0.5)\n");
    printf("  -y <double>   Y center (default=0.0)\n");
    printf("  -s <double>   Scale (default=4.0)\n");
    printf("  -W <int>      Image width (default=1000)\n");
    printf("  -H <int>      Image height (default=1000)\n");
    printf("  -m <int>      Max iterations (default=1000)\n");
    printf("  --movie       Generate a movie instead of an image\n");
}