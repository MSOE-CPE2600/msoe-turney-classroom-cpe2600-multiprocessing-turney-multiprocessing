//
//  mandel.c
//  Name: Alion Bujku
//  Enhanced to support multiprocessing (-p) and multithreading (-t) simultaneously.
//

// Standard library headers for general operations
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/wait.h>
#include <pthread.h> // For multithreading
#include <semaphore.h> // For inter-process synchronization
#include <fcntl.h> // For file control options
#include <math.h> // For mathematical functions
#include <getopt.h> // For parsing command-line options
#include "jpegrw.h" // Custom library for handling JPEG files

// Struct to encapsulate data passed to threads
typedef struct {
    imgRawImage *img; // Pointer to the image being generated
    double xmin, xmax, ymin, ymax; // Coordinates for the Mandelbrot calculation
    int max; // Maximum number of iterations
    int y_start, y_end; // Range of rows to compute
} ThreadData;

// Function prototypes
int iterations_at_point(double x, double y, int max); // Calculate iterations for a point
int iteration_to_color(int i, int max); // Map iteration count to a color
void compute_image_multithreaded(imgRawImage *img, double xmin, double xmax, double ymin, double ymax, int max, int num_threads);
void *thread_compute(void *args); // Worker thread function for Mandelbrot computation
void generate_movie(double xcenter, double ycenter, double xscale, int image_width, int image_height, int max, int num_frames, int num_processes, int num_threads);
void show_help(); // Display help message

int main(int argc, char *argv[]) {
    // Default parameters for Mandelbrot image generation
    double xcenter = -0.5, ycenter = 0.0; // Default center of the fractal
    double xscale = 4.0, yscale = 0.0; // Scale of the fractal (yscale computed later)
    int image_width = 1000, image_height = 1000; // Image resolution
    int max = 1000; // Maximum iterations for Mandelbrot calculation
    int num_threads = 1, num_processes = 1, num_frames = 50; // Defaults for parallelism and movie frames
    int is_movie_mode = 0; // Flag to indicate movie generation mode
    const char *outfile = NULL; // Output file name for single image mode

    // Command-line options definition
    struct option long_options[] = {
        {"movie", no_argument, &is_movie_mode, 1}, // Enable movie mode
        {"x", required_argument, 0, 'x'}, // X center coordinate
        {"y", required_argument, 0, 'y'}, // Y center coordinate
        {"scale", required_argument, 0, 's'}, // Fractal scale
        {"width", required_argument, 0, 'W'}, // Image width
        {"height", required_argument, 0, 'H'}, // Image height
        {"max", required_argument, 0, 'm'}, // Max iterations
        {"frames", required_argument, 0, 'f'}, // Number of frames for movie
        {"processes", required_argument, 0, 'p'}, // Number of processes
        {"threads", required_argument, 0, 't'}, // Number of threads
        {"output", required_argument, 0, 'o'}, // Output file name
        {"help", no_argument, 0, 'h'}, // Show help
        {0, 0, 0, 0} // End of options
    };

    // Parse command-line arguments
    int opt, option_index = 0;
    while ((opt = getopt_long(argc, argv, "x:y:s:W:H:m:o:f:p:t:h", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'x': xcenter = atof(optarg); break; // Set X center
            case 'y': ycenter = atof(optarg); break; // Set Y center
            case 's': xscale = atof(optarg); break; // Set scale
            case 'W': image_width = atoi(optarg); break; // Set image width
            case 'H': image_height = atoi(optarg); break; // Set image height
            case 'm': max = atoi(optarg); break; // Set max iterations
            case 'o': outfile = optarg; break; // Set output file name
            case 'f': num_frames = atoi(optarg); break; // Set number of movie frames
            case 'p': num_processes = atoi(optarg); break; // Set number of processes
            case 't': num_threads = atoi(optarg); break; // Set number of threads
            case 'h': show_help(); return 0; // Show help message
            default: show_help(); return 1; // Invalid input
        }
    }

    // Compute yscale to match aspect ratio
    yscale = xscale / image_width * image_height;

    // Handle movie generation mode
    if (is_movie_mode) {
        generate_movie(xcenter, ycenter, xscale, image_width, image_height, max, num_frames, num_processes, num_threads);
        return 0;
    }

    // Validate single image mode (output file is mandatory)
    if (outfile == NULL) {
        fprintf(stderr, "Error: Output file (-o) must be specified for single image mode.\n");
        return 1;
    }

    // Initialize image structure
    imgRawImage *img = initRawImage(image_width, image_height);

    // Compute Mandelbrot set using multithreading
    compute_image_multithreaded(img, xcenter - xscale / 2, xcenter + xscale / 2,
                                ycenter - yscale / 2, ycenter + yscale / 2, max, num_threads);

    // Save the final image
    if (storeJpegImageFile(img, outfile) != 0) {
        fprintf(stderr, "Error: Failed to save image to %s\n", outfile);
    }

    // Free allocated resources
    freeRawImage(img);

    return 0;
}

// Generate a Mandelbrot movie using multiprocessing and multithreading
void generate_movie(double xcenter, double ycenter, double xscale, int image_width, int image_height, int max, int num_frames, int num_processes, int num_threads) {
    sem_unlink("/mandel_sem"); // Unlink existing semaphore
    sem_t *sem = sem_open("/mandel_sem", O_CREAT | O_EXCL, 0644, num_processes); // Initialize semaphore
    if (sem == SEM_FAILED) {
        perror("sem_open failed");
        exit(1);
    }

    double zoom_factor = 0.90; // Zoom factor for each frame

    for (int frame = 0; frame < num_frames; frame++) {
        sem_wait(sem); // Wait for available process slot

        if (fork() == 0) { // Child process
            char filename[256];
            snprintf(filename, sizeof(filename), "frame%03d.jpg", frame);

            double scale = xscale * pow(zoom_factor, frame);
            double xmin = xcenter - scale / 2;
            double xmax = xcenter + scale / 2;
            double ymin = ycenter - scale / 2;
            double ymax = ycenter + scale / 2;

            imgRawImage *img = initRawImage(image_width, image_height);
            compute_image_multithreaded(img, xmin, xmax, ymin, ymax, max, num_threads);
            storeJpegImageFile(img, filename);
            freeRawImage(img);

            sem_post(sem); // Signal semaphore
            exit(0); // End child process
        }
    }

    while (wait(NULL) > 0); // Wait for all child processes

    sem_close(sem); // Close semaphore
    sem_unlink("/mandel_sem"); // Unlink semaphore
}

// Mandelbrot computation using multithreading
void compute_image_multithreaded(imgRawImage *img, double xmin, double xmax, double ymin, double ymax, int max, int num_threads) {
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];
    int rows_per_thread = img->height / num_threads;

    for (int i = 0; i < num_threads; i++) {
        thread_data[i] = (ThreadData){
            .img = img,
            .xmin = xmin,
            .xmax = xmax,
            .ymin = ymin,
            .ymax = ymax,
            .max = max,
            .y_start = i * rows_per_thread,
            .y_end = (i == num_threads - 1) ? img->height : (i + 1) * rows_per_thread
        };
        pthread_create(&threads[i], NULL, thread_compute, &thread_data[i]);
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
}

// Worker thread function for Mandelbrot computation
void *thread_compute(void *args) {
    ThreadData *data = (ThreadData *)args;

    for (int y = data->y_start; y < data->y_end; y++) {
        for (int x = 0; x < data->img->width; x++) {
            double cx = data->xmin + x * (data->xmax - data->xmin) / (data->img->width - 1);
            double cy = data->ymin + y * (data->ymax - data->ymin) / (data->img->height - 1);
            int iterations = iterations_at_point(cx, cy, data->max);
            int color = iteration_to_color(iterations, data->max);
            setPixelCOLOR(data->img, x, y, color);
        }
    }

    return NULL;
}

// Calculate Mandelbrot iterations for a given point
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

// Map iteration count to an RGB color
int iteration_to_color(int i, int max) {
    if (i == max) return 0x000000; // Black for points inside the set
    int r = (i * 9) % 256;
    int g = (i * 15) % 256;
    int b = (i * 25) % 256;
    return (r << 16) | (g << 8) | b;
}

// Display help message
void show_help() {
    printf("Usage: mandel [options]\n");
    printf("Options:\n");
    printf("  --movie       Generate a movie instead of a single image\n");
    printf("  -x <double>   X center (default=-0.5)\n");
    printf("  -y <double>   Y center (default=0.0)\n");
    printf("  -s <double>   Scale (default=4.0)\n");
    printf("  -W <int>      Image width (default=1000)\n");
    printf("  -H <int>      Image height (default=1000)\n");
    printf("  -m <int>      Max iterations (default=1000)\n");
    printf("  -o <string>   Output file name (required for single image mode)\n");
    printf("  -f <int>      Number of frames for movie (default=50)\n");
    printf("  -p <int>      Number of processes (default=1)\n");
    printf("  -t <int>      Number of threads per process (default=1)\n");
    printf("  -h            Show this help message\n");
}