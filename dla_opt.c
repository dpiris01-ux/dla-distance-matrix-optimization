#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

unsigned int xorshift_state;

// Función para inicializar el generador con una semilla

void xorshift_seed(unsigned int seed) {
    if (seed == 0) {
        seed = time(NULL); // Evitar semilla 0 que puede causar problemas
    }
    xorshift_state = seed;
}

// Generador Xorshift que devuelve un entero de 32 bits

unsigned int xorshift32() {
    unsigned int x = xorshift_state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    xorshift_state = x;
    return x;
}

// Función que devuelve un double aleatorio de alta calidad en el rango [0.0, 1.0)

long double rand_double() {
    // Se divide por 2^32
    return xorshift32() / 4294967296.0;
}


// Parámetros de la simulación

#define L 8192     // Tamaño de la caja
#define T 100000    // Número de partículas a agregar
#define N 1        // Número de muestras 
#define layers 100   // Capas para análisis intermedio
#define PARTICLES_PER_LAYER (T/layers)

// Parámetros de la partícula y la caminata

#define diameter 2.0

// Radios de lanzamiento y eliminación

#define rmax 5

// Estructuras de datos de la optimización

#define L_PARAM 4 // l
#define CELL_SIZE L_PARAM // Tamaño de la celda, en este caso l
#define GRID_DIM (L/CELL_SIZE) // Dimensión de grilla de tamaño lxl
#define BLOCK_SIZE (L_PARAM*L_PARAM) // Tamaño fijo de l^2 para cada bloque de celdas

#define MAX_GRID_CELLS (GRID_DIM*GRID_DIM)

// Matriz W: Mapea coordenadas (i,j) a un índice secuencial k. 0 = celda no usada.

int W_grid[GRID_DIM][GRID_DIM];
int next_k_value; // Siguiente k a asignar

// Matriz Nk: Almacena el número de partículas en la celda k.

int Nk_counts[MAX_GRID_CELLS + 1];

// Matriz F: El array 1D que contiene los índices de las partículas en bloques contiguos.

int F_indices[MAX_GRID_CELLS*BLOCK_SIZE];

// Optimización del análisis de interfaz

#define NUM_REGIONES 2*(36) // Por ejemplo, 36 regiones de 10 grados cada una.
#define PARTICLES_TO_ANALYZE_PER_LAYER (PARTICLES_PER_LAYER / 10)
#define LAST_PARTICLES_TO_SAVE 100

void add_particle_to_grid(int p_index);
int check_collision(long double walker_x, long double walker_y);
int find_particle_at(long double search_x, long double search_y);

long double max_radios_por_region[NUM_REGIONES];
long double active_zone_radii[PARTICLES_TO_ANALYZE_PER_LAYER];

// Variables globales

int cx, cy, stuck, walking, rt_count, particle_count, i, j, t, m,
particle_hit, PARTICLES[layers], num_layers, step_index, step_index_final, r2, x_idx, y_idx,
d, particle_found_flags[T+1];

double particles_list[T + 1][2], max_r, walker_x, walker_y, launch_radius, theta,
move_theta, current_r, escape_dist_sq, rt[T], angular_step, theta_r, search_x, search_y,
px, py, dist_sq, search_step, current_radius, scan_x, scan_y, collision_dist_sq,
total_rg_per_layer[layers];

// Interfaz

long double skewnessA[layers], kurtosisA[layers], varA_vector[layers], k3B_vector[layers], k4B_vector[layers], varB_vector[layers],
varA, skA, ktA, varB, skB, ktB, varC, skC, ktC, k3A, k4A, k3B, k4B, k3C, k4C, sqr_varB, skewnessC, kurtosisC, skewnessB, kurtosisB;

// Zona activa

long double skewnessA_active[layers], kurtosisA_active[layers], varA_vector_active[layers], k3B_vector_active[layers],
k4B_vector_active[layers], varB_vector_active[layers], varA_active, skA_active, ktA_active, varB_active, skB_active, ktB_active,
k3A_active, k4A_active, skewness_activeB, kurtosis_activeB;

// Radio máximo

long double skewnessA_rmax[layers], kurtosisA_rmax[layers], varA_vector_rmax[layers], k3B_vector_rmax[layers], k4B_vector_rmax[layers],
varB_vector_rmax[layers], varA_rmax, skA_rmax, ktA_rmax, varB_rmax, skB_rmax, ktB_rmax, k3A_rmax, k4A_rmax, k4B_rmax, k3B_rmax,
skewness_rmaxB, kurtosis_rmaxB;

// Todas las partículas

double radii_list[T + 1];
long double varA_vector_all[layers], skewnessA_all[layers], kurtosisA_all[layers], varB_vector_all[layers], k3B_vector_all[layers],
k4B_vector_all[layers], M1_all[layers], M2_all[layers], M3_all[layers], M4_all[layers];

long double M1[layers], M2[layers], M3[layers], M4[layers], varC, ri, mr1, mr2, mr3, mr4, total_results[60][2], theta_radial, skewnessC,
kurtosisC ,k3C, k4C, M1_active[layers], M2_active[layers], M3_active[layers], M4_active[layers], num_rt_particles[L],
M1_rmax[layers], M2_rmax[layers], M3_rmax[layers], M4_rmax[layers];

// Variables para la distribución del radio media, la rugosidad y el radio de giro

long double M1_mean_dist[layers], M2_mean_dist[layers], M3_mean_dist[layers], M4_mean_dist[layers], M1_rough_dist[layers],
M2_rough_dist[layers], M3_rough_dist[layers], M4_rough_dist[layers], M1_giro_dist[layers], M2_giro_dist[layers],
M3_giro_dist[layers], M4_giro_dist[layers];

long double M1_rmax_dist[layers], M2_rmax_dist[layers], M3_rmax_dist[layers], M4_rmax_dist[layers];


// Distancia máxima que la grid Ω calculará con precisión.

#define D_MAX 90

// Paso mínimo del caminante justo antes de una colisión.

#define L_MIN 1.0

// Condición para el chequeo de colisión detallado (2*r_p + L_min + 1)

#define COLLISION_THRESHOLD 4


// Red Y (On-lattice cluster grid)

unsigned short int Y_grid[L][L];

// Red Ω (On-lattice distance grid)

unsigned short int Omega_grid[L][L];

// Red Ψ (Vicinity grid)

unsigned short int Vicinity_grid[2 * D_MAX + 1][2 * D_MAX + 1];

// Función para mostrar barra de progreso

void mostrar_barra_progreso(int actual, int total, const char* prefijo, time_t tiempo_inicio) {

    int porcentaje = (actual * 100)/total;

    // Calcular tiempo transcurrido

    time_t tiempo_actual = time(NULL);
    int segundos = (int)difftime(tiempo_actual, tiempo_inicio);
    int minutos = segundos/60;
    int segundos_restantes = segundos % 60;

    // Limpiar la línea anterior y sobrescribir

    printf("\r\033[K");

    // Imprimir barra de progreso con tiempo en formato mm:ss

    printf("%s: [", prefijo);
    for (int j = 0; j < 60; j++) {
        if (j < porcentaje/2) printf("#");
        else printf(" ");
    }
    printf("] %3d%% | Tiempo: %02d:%02d", porcentaje, minutos, segundos_restantes);
    fflush(stdout);
}

// Función para inicializar todas las nuevas redes al comienzo de una simulación.

void inicializar_grids() {

    // Calcular la Vicinity Grid (Ψ) una sola vez
    // Esta grid es una plantilla de distancias a su propio centro.

    for (int i = 0; i < 2 * D_MAX + 1; i++) {

        for (int j = 0; j < 2 * D_MAX + 1; j++) {

            // Calcular distancia desde (i, j) al centro (D_MAX, D_MAX)

            long double dist = sqrt((i - D_MAX)*(i - D_MAX) + (j - D_MAX)*(j - D_MAX));

            Vicinity_grid[i][j] = (int)round(dist);
        }
    }

    // Inicializar Grids Y y Ω
    // La red de distancias (Ω) se llena con D_MAX, y la red de clúster (Y) con 0 (vacío).

    for (int i = 0; i < L; i++) {

        for (int j = 0; j < L; j++) {

            Omega_grid[i][j] = D_MAX;
            Y_grid[i][j] = 0;
        }
    }
}

// Coloca la partícula semilla y actualiza las redes por primera vez.

void colocar_particula_semilla() {

    // La semilla ahora está en el origen (0,0) del espacio de simulación.

    particles_list[0][0] = 0.0;
    particles_list[0][1] = 0.0;
    particle_count = 1;
    max_r = 1.5*diameter;
    
    // El centro de la GRID es (L/2, L/2)

    int seed_grid_x = L/2;
    int seed_grid_y = L/2;

    add_particle_to_grid(0); 
    
    Y_grid[seed_grid_x][seed_grid_y] = 1; // El índice 1 en el centro de la red

    // Actualiza Omega_grid alrededor del centro de la red

    for (int i = -D_MAX; i <= D_MAX; i++) {

        for (int j = -D_MAX; j <= D_MAX; j++) {

            int grid_x = seed_grid_x + i;
            int grid_y = seed_grid_y + j;

            if (grid_x >= 0 && grid_x < L && grid_y >= 0 && grid_y < L) Omega_grid[grid_x][grid_y] = Vicinity_grid[i + D_MAX][j + D_MAX];
        }
    }
}

// Actualiza las redes Y y Ω después de que una nueva partícula se ha adherido.

void actualizar_grids_por_nueva_particula(int p_index) {

    int px_grid = (int)round((L / 2.0) + particles_list[p_index][0]);
    int py_grid = (int)round((L / 2.0) + particles_list[p_index][1]);

    // Usa las nuevas variables para los índices de la red

    if (px_grid < 0 || px_grid >= L || py_grid < 0 || py_grid >= L) return;

    Y_grid[px_grid][py_grid] = p_index + 1;

    for (int i = -D_MAX; i <= D_MAX; i++) {

        int gx = px_grid + i;
        if (gx < 0 || gx >= L) continue;

        for (int j = -D_MAX; j <= D_MAX; j++) {

            int gy = py_grid + j;
            if (gy < 0 || gy >= L) continue;
            
            int d = Vicinity_grid[i + D_MAX][j + D_MAX];

            if (d < Omega_grid[gx][gy]) Omega_grid[gx][gy] = d;
        }
    }
}

// Funciones para la optimización de la red

void initialize_coarse_grid() {

    // Inicializar W, donde 0 significa "sin asignar"

    for (int i = 0; i < GRID_DIM; i++) {

        for (int j = 0; j < GRID_DIM; j++) {

            W_grid[i][j] = 0;
        }
    }

    next_k_value = 1; // El primer índice secuencial será 1
    
    // Inicializar contadores de partículas para todas las posibles celdas k
    // (usar MAX_GRID_CELLS+1 porque k empieza en 1)

    for (int i = 0; i <= MAX_GRID_CELLS; i++) Nk_counts[i] = 0;

    // Inicializar F_indices a -1 para depuración

    for (int i = 0; i < MAX_GRID_CELLS * BLOCK_SIZE; i++) F_indices[i] = -1;

}

void add_particle_to_grid(int p_index) {

    double px = particles_list[p_index][0];
    double py = particles_list[p_index][1];

    double world_to_grid_offset = (long double)(L/2);

    int cell_x = (int)round((px + world_to_grid_offset)/CELL_SIZE);
    int cell_y = (int)round((py + world_to_grid_offset)/CELL_SIZE);
    if (cell_x < 0 || cell_x >= GRID_DIM || cell_y < 0 || cell_y >= GRID_DIM) return;

    // Asignar el índice secuencial 'k' de la celda

    int k = W_grid[cell_x][cell_y];

    if (k == 0) { // Si la celda se ocupa por primera vez

        k = next_k_value;
        W_grid[cell_x][cell_y] = k;
        next_k_value++;

    }

    // Obtener el número actual de partículas en la celda k

    int count = Nk_counts[k];

    // Calcular la posición en el array F para añadir el nuevo índice

    int start_of_block = BLOCK_SIZE * (k - 1);
    int position_to_add = start_of_block + count;

    F_indices[position_to_add] = p_index;

    // Incrementar el contador de partículas para la celda k

    Nk_counts[k]++;
}

int check_collision(long double walker_x, long double walker_y) {

    int cell_x = (int)round((walker_x + (L/2.0))/CELL_SIZE);
    int cell_y = (int)round((walker_y + (L/2.0))/CELL_SIZE);

    // Revisar la celda del walker y sus 8 vecinas

    for (int dx = -1; dx <= 1; dx++) {

        for (int dy = -1; dy <= 1; dy++) {

            int check_cell_x = cell_x + dx;
            int check_cell_y = cell_y + dy;

            if (check_cell_x >= 0 && check_cell_x < GRID_DIM && check_cell_y >= 0 && check_cell_y < GRID_DIM) {
                
                int k = W_grid[check_cell_x][check_cell_y];
                if (k == 0) continue; // Celda vacía, no hay nada que revisar

                int num_particles_in_cell = Nk_counts[k];
                int start_of_block = BLOCK_SIZE * (k - 1);

                // Recorrer el bloque de memoria contiguo en F 

                for (int i = num_particles_in_cell - 1; i >= 0; i--) {

                    int p_index = F_indices[start_of_block + i];

                    double cluster_p_x = particles_list[p_index][0];
                    double cluster_p_y = particles_list[p_index][1];

                    double dist_x = walker_x - cluster_p_x;
                    double dist_y = walker_y - cluster_p_y;
                    double distance_sq = dist_x*dist_x + dist_y*dist_y;

                    if (distance_sq < collision_dist_sq) return 1;
                }
            }
        }
    }

    return 0;
}

// Función optimización para búsqueda 

int find_particle_at(long double search_x, long double search_y) {

    int cell_x = (int)round((search_x + (L/2.0))/CELL_SIZE);
    int cell_y = (int)round((search_y + (L/2.0))/CELL_SIZE);

    // Revisar la celda del punto de búsqueda y sus 8 vecinas

    for (int dx = -1; dx <= 1; dx++) {

        for (int dy = -1; dy <= 1; dy++) {

            int check_cell_x = cell_x + dx;
            int check_cell_y = cell_y + dy;

            if (check_cell_x >= 0 && check_cell_x < GRID_DIM && check_cell_y >= 0 && check_cell_y < GRID_DIM) {

                int k = W_grid[check_cell_x][check_cell_y];
                if (k == 0) continue; 

                int num_particles_in_cell = Nk_counts[k];
                int start_of_block = BLOCK_SIZE * (k - 1);

                for (int i = 0; i < num_particles_in_cell; i++) {

                    int p_index = F_indices[start_of_block + i];
                    double p_x = particles_list[p_index][0];
                    double p_y = particles_list[p_index][1];

                    long double dist_sq = (search_x - p_x)*(search_x - p_x) + (search_y - p_y)*(search_y - p_y);

                    if (dist_sq < collision_dist_sq) return p_index;
                }
            }
        }
    }

    return -1; // No se encontró ninguna partícula cerca
}

// Función para guardar partícula

void save_particles(const char* filename) {
    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("No se pudo abrir el archivo para guardar partículas");
        return;
    }
    for (int i = 0; i < particle_count; i++) {
        fprintf(file, "%.6f %.6f\n", particles_list[i][0], particles_list[i][1]);
    }
    fclose(file);
}

// Función para guardar la distribución de radios

void save_radii_distribution(const char* filename, int count) {

    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("No se pudo abrir el archivo para guardar los radios");
        return;
    }
    for (int i = 0; i < count; i++) {
        fprintf(file, "%.6f\n", rt[i]);
    }
    fclose(file);
}

// Función para calcular el radio de giro al cuadrado

double get_rg_sq(int current_particle_count) {

    double sum_x = 0.0, sum_y = 0.0, sum_x_sq = 0.0, sum_y_sq = 0.0;

    for (int i = 0; i < current_particle_count; i++) {

        double px = particles_list[i][0];
        double py = particles_list[i][1];
        
        sum_x += px;
        sum_y += py;
        sum_x_sq += px*px;
        sum_y_sq += py*py;
    }

    double N_p = (double)current_particle_count;
    double mean_x = sum_x/N_p;
    double mean_y = sum_y/N_p;
    double mean_x_sq = sum_x_sq/N_p;
    double mean_y_sq = sum_y_sq/N_p;

    double rg_sq = (mean_x_sq + mean_y_sq) - (mean_x*mean_x + mean_y*mean_y);
    
    return (rg_sq > 0) ? rg_sq : 0.0;
}


// Calcula un "pseudo-ángulo" para el vector (dx, dy).

// El valor de retorno está en el rango [0, 8) y aumenta a medida que el ángulo rota

double pseudo_angulo(double dx, double dy) {

    double ax = fabs(dx);
    double ay = fabs(dy);
    double p = (ax > ay) ? ay / ax : ax / ay;

    if (dx >= 0 && dy >= 0) { // Cuadrante 1
        return (ax > ay) ? p : 2.0 - p; // Octantes 0 y 1
    } else if (dx < 0 && dy >= 0) { // Cuadrante 2
        return (ax > ay) ? 4.0 - p : 2.0 + p; // Octantes 3 y 2
    } else if (dx < 0 && dy < 0) { // Cuadrante 3
        return (ax > ay) ? 4.0 + p : 6.0 - p; // Octantes 4 y 5
    } else { // Cuadrante 4 (dx >= 0 && dy < 0)
        return (ax > ay) ? 8.0 - p : 6.0 + p; // Octantes 7 y 6
    }
}

// Implementación del colisión

// Resuelve la ecuación cuadrática a*x^2 + b*x + c = 0
// Devuelve la raíz más pequeña y positiva, o -1.0 si no existe tal raíz.

double solve_quadratic(double a, double b, double c) {

    double discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0) {
        // No hay soluciones reales, no hay colisión en esta trayectoria.
        return -1.0;
    }

    double sqrt_discriminant = sqrt(discriminant);
    double root1 = (-b + sqrt_discriminant)/(2.0 * a);
    double root2 = (-b - sqrt_discriminant)/(2.0 * a);

    // Nos interesa la colisión más cercana que ocurra en el futuro (raíz positiva más pequeña).
    // Usamos una pequeña tolerancia (epsilon) para evitar problemas con el punto flotante cerca de cero.

    double epsilon = 1e-9;

    if (root1 > epsilon && root2 > epsilon) {
        return fmin(root1, root2);
    } else if (root1 > epsilon) {
        return root1;
    } else if (root2 > epsilon) {
        return root2;
    }

    // Ambas raíces son negativas o cero, la colisión está en el pasado.

    return -1.0;
}

// Encuentra la distancia de colisión más cercana (L_hit) para un caminante.
// Devuelve la L_hit más pequeña, o -1.0 si no hay colisión dentro del paso L_max.

double find_collision_distance(double walker_x, double walker_y, double move_theta, double L_max) {
    
    double min_L_hit = L_max + 1.0; // Inicializar con un valor mayor que L_max
    
    // Convertir la posición del caminante a índices de la grilla gruesa

    int cell_x = (int)round((walker_x + (L/2.0))/CELL_SIZE);
    int cell_y = (int)round((walker_y + (L/2.0))/CELL_SIZE);
    
    double cos_theta = cos(move_theta);
    double sin_theta = sin(move_theta);

    // Revisar la celda del walker y sus 8 vecinas

    for (int dx = -1; dx <= 1; dx++) {

        for (int dy = -1; dy <= 1; dy++) {

            int check_cell_x = cell_x + dx;
            int check_cell_y = cell_y + dy;

            if (check_cell_x >= 0 && check_cell_x < GRID_DIM && check_cell_y >= 0 && check_cell_y < GRID_DIM) {

                int k = W_grid[check_cell_x][check_cell_y];
                if (k == 0) continue; // Celda vacía

                int num_particles_in_cell = Nk_counts[k];
                int start_of_block = BLOCK_SIZE*(k - 1);

                // Iterar sobre cada partícula en la celda candidata

                for (int i = num_particles_in_cell - 1; i >= 0; i--) {

                    int p_index = F_indices[start_of_block + i];
                    double p_x = particles_list[p_index][0];
                    double p_y = particles_list[p_index][1];

                    double max_safe_dist = L_max + diameter;

                    if (fabs(walker_x - p_x) > max_safe_dist || fabs(walker_y - p_y) > max_safe_dist) {
                        continue; // Saltar a la siguiente partícula
                    }

                    // Calcular coeficientes A, B, C 

                    double B = 2.0*(cos_theta*(walker_x - p_x) + sin_theta*(walker_y - p_y));

                    // Collision_dist_sq es diameter*diameter

                    double C = (walker_x - p_x)*(walker_x - p_x) + (walker_y - p_y)*(walker_y - p_y) - collision_dist_sq;

                    double L_hit = solve_quadratic(1.0, B, C);

                    // Si se encontró una colisión válida (0 < L_hit <= L_max) y es más cercana que la anterior, se guarda

                    if (L_hit > 1e-9 && L_hit <= L_max) {

                        if (L_hit < min_L_hit) min_L_hit = L_hit;
                    }
                }
            }
        }
    }

    if (min_L_hit <= L_max) {
        return min_L_hit;
    } else {
        return -1.0; // No se encontró colisión
    }
}

// Inicio de la simulación

int main(int argc, char *argv[]) {

    int job_id = 0; // ID por defecto

    if (argc > 1) {
        job_id = atoi(argv[1]); // Capturar el ID desde la línea de comandos
    }

    // Garantizar que incluso si dos trabajos empiezan en el mismo segundo, la semilla será diferente.

    xorshift_seed(time(NULL) ^ getpid() ^ job_id); 
    
    time_t tiempo_inicio = time(NULL);
    
    // Calcular una sola vez la distancia de colisión al cuadrado (más rápido que hacer sqrt)

    collision_dist_sq = diameter * diameter;

    for (int layer_ini = 0; layer_ini < layers; layer_ini++){

            M1_mean_dist[layer_ini] = 0;
            M2_mean_dist[layer_ini] = 0; 
            M3_mean_dist[layer_ini] = 0;
            M4_mean_dist[layer_ini] = 0; 

            M1_rough_dist[layer_ini] = 0;
            M2_rough_dist[layer_ini] = 0;
            M3_rough_dist[layer_ini] = 0;
            M4_rough_dist[layer_ini] = 0;

            M1_giro_dist[layer_ini] = 0;
            M2_giro_dist[layer_ini] = 0;
            M3_giro_dist[layer_ini] = 0;
            M4_giro_dist[layer_ini] = 0;

        }

    cx = 0;
    cy = 0;

    for (m = 0; m < N; m++) {

        mostrar_barra_progreso(m + 1, N, "Procesando Muestra", tiempo_inicio);


        initialize_coarse_grid(); // Grid para check_collision
        inicializar_grids();      // Las nuevas redes
        colocar_particula_semilla(); // Coloca la partícula 0 y actualiza las redes
        

        for (int i = 0; i < NUM_REGIONES; i++) max_radios_por_region[i] = 0.0;
        int active_zone_radii_count = 0;

        // Bucle único para todas las partículas

        for (t = 1; t < T; t++) {

            stuck = 0;

            while (!stuck) {

                // Lanza el caminante relativo al origen (0,0)

                long double launch_radius = max_r + rmax;
                long double theta = rand_double()*2.0*M_PI;

                walker_x = launch_radius*cos(theta);
                walker_y = launch_radius*sin(theta);

                walking = 1;

                while (walking) {

                    long double current_step;

                    // Obtener la posición del caminante en la red de distancias

                    int walker_ix = (L/2) + (int)round(walker_x);
                    int walker_iy = (L/2) + (int)round(walker_y);


                    // Consultar la distancia estimada al clúster desde la red omega

                    int d_wc = D_MAX; // Valor por defecto si está fuera de la red
                    if (walker_ix >= 0 && walker_ix < L && walker_iy >= 0 && walker_iy < L) d_wc = Omega_grid[walker_ix][walker_iy];

                    // Decidir el tamaño del paso según la distancia

                    if (d_wc <= COLLISION_THRESHOLD) {

                        // Está muy cerca, posible colisión. Usar paso mínimo

                        current_step = L_MIN;

                    } else if (d_wc < D_MAX) {

                        // Está a una distancia intermedia segura. Dar un salto grande

                        current_step = d_wc - COLLISION_THRESHOLD;

                    } else {

                        // Está muy lejos (d_wc == D_MAX). Dar el salto más grande posible

                        current_step = D_MAX - COLLISION_THRESHOLD;

                    }

                    // Mover la partícula

                    long double move_theta = rand_double()*2.0*M_PI;

                    if (current_step == L_MIN) {

                        // Usamos la nueva función para predecir la colisión

                        double L_hit = find_collision_distance(walker_x, walker_y, move_theta, L_MIN);

                        if (L_hit > 0) { // Se encontró una colisión válida

                            // Mover el caminante la distancia exacta para la colisión

                            walker_x += L_hit * cos(move_theta);
                            walker_y += L_hit * sin(move_theta);

                            particles_list[particle_count][0] = walker_x;
                            particles_list[particle_count][1] = walker_y;
                            
                            add_particle_to_grid(particle_count);
                            actualizar_grids_por_nueva_particula(particle_count);
                            
                            particle_count++;

                            double current_r = sqrt(walker_x*walker_x + walker_y*walker_y); // Asumiendo cx,cy = 0

                           if (current_r > max_r) max_r = current_r;

                            radii_list[particle_count] = current_r;
                            
                            // Calcular el ángulo de la nueva partícula y el vector desde el centro

                            double dx = walker_x - cx;
                            double dy = walker_y - cy;

                            // Obtiene el pseudo-ángulo (rango 0-8)

                            double valor_pseudo = pseudo_angulo(dx, dy);

                            // Escala el resultado del pseudo-ángulo al número de regiones

                            // El rango de valor_pseudo es [0, 8), por lo que dividimos por 8.0

                            int region_idx = (int)(valor_pseudo*(NUM_REGIONES/8.0));

                            // Asegurarse de que el índice esté dentro de los límites por si acas
                            if (region_idx >= NUM_REGIONES) region_idx = NUM_REGIONES - 1;

                            // Actualiza el radio máximo para ese sector

                            if (current_r > max_radios_por_region[region_idx]) max_radios_por_region[region_idx] = current_r;

                            // Si estamos en la segunda mitad de las partículas de esta capa

                            int particles_in_current_layer = t%PARTICLES_PER_LAYER;
                        
                            if (particles_in_current_layer >= PARTICLES_PER_LAYER - LAST_PARTICLES_TO_SAVE) {

                                // Guardamos el radio de esta partícula para analizarlo después

                                if (active_zone_radii_count < PARTICLES_TO_ANALYZE_PER_LAYER) {

                                    active_zone_radii[active_zone_radii_count] = current_r;
                                    active_zone_radii_count++;
                                }
                            }
                            
                            stuck = 1;
                            walking = 0;

                        } else { // No hubo colisión, dar el paso completo

                            walker_x += current_step*cos(move_theta);
                            walker_y += current_step*sin(move_theta);

                        }

                    } else { // El paso era grande (lejos del clúster), no es necesario chequear colisión

                        walker_x += current_step*cos(move_theta);
                        walker_y += current_step*sin(move_theta);
                    }

                    double escape_dist_sq = (walker_x - cx)*(walker_x - cx) + (walker_y - cy)*(walker_y - cy);

                    // Radio de muerte 5*Rmax

                    if (escape_dist_sq > (5*max_r)*(5*max_r)) walking = 0;

                }

            } // Fin while walking 1

            if ((particle_count % PARTICLES_PER_LAYER) == 0 && particle_count > 0) {
                
                int layer_idx = (particle_count/PARTICLES_PER_LAYER) - 1; 

                // Análisis de la interfaz

                rt_count = 0;

                for (int flag = 0; flag < T+1; flag++) particle_found_flags[flag] = 0;
                
                angular_step = diameter/max_r;

                for (theta_r = 0; theta_r < 2.0*M_PI; theta_r += angular_step) {

                    int region_idx = (int)(theta_r / (2.0*M_PI)*NUM_REGIONES);
                    double region_max_r = max_radios_por_region[region_idx];
                    
                    search_step = diameter/2.0;

                    for (current_radius = region_max_r + diameter; current_radius > 0; current_radius -= search_step) {
                       scan_x = cx + current_radius*cos(theta_r);
                        scan_y = cy + current_radius*sin(theta_r);

                        int found_p_index = find_particle_at(scan_x, scan_y);

                        if (found_p_index != -1) {

                            if (!particle_found_flags[found_p_index]) {

                                particle_found_flags[found_p_index] = 1;

                                double px = particles_list[found_p_index][0];
                                double py = particles_list[found_p_index][1];

                                rt[rt_count] = sqrt((px - cx)*(px - cx) + (py - cy)*(py - cy));
                                rt_count++;
                            }

                            break;
                        }
                    }

                }


                ri = 0, mr1 = 0, mr2 = 0, mr3 = 0, mr4 = 0;

                long double mr1_active = 0, mr2_active = 0, mr3_active = 0, mr4_active = 0;

                long double rmax1= 0, rmax2 = 0, rmax3 = 0, rmax4 = 0;

                for (int mean_idx = 0; mean_idx < rt_count; mean_idx++){

                    ri = (long double)rt[mean_idx];
                    mr1 += ri;
                    mr2 += ri*ri;
                    mr3 += ri*ri*ri;
                    mr4 += ri*ri*ri*ri;

                }

                mr1 /= (long double)rt_count;
                mr2 /= (long double)rt_count;
                mr3 /= (long double)rt_count;
                mr4 /= (long double)rt_count;

                M1[layer_idx] += mr1;
                M2[layer_idx] += mr2;
                M3[layer_idx] += mr3;
                M4[layer_idx] += mr4;

                varA = mr2 - (mr1*mr1);
                k3A = mr3 - 3*mr1*mr2 + 2*(mr1*mr1*mr1);
                k4A = mr4 - 4*mr1*mr3 - 3*(mr2*mr2) + 12*(mr1*mr1)*mr2 - 6*(mr1*mr1*mr1*mr1);

                varA_vector[layer_idx] += varA;
                skewnessA[layer_idx] += k3A/(varA*sqrt(varA));
                kurtosisA[layer_idx] += k4A/(varA*varA);

                varB_vector[layer_idx] += mr2 - (mr1*mr1);
                k3B_vector[layer_idx] += mr3 - 3*mr1*mr2 + 2*(mr1*mr1*mr1);
                k4B_vector[layer_idx] += mr4 - 4*mr1*mr3 - 3*(mr2*mr2) + 12*(mr1*mr1)*mr2 - 6*(mr1*mr1*mr1*mr1);

                num_rt_particles[layer_idx] += (long double)rt_count;

                // Utilizando el método de Paul Meakin

                for (int mean_idx = 0; mean_idx < active_zone_radii_count; mean_idx++) {

                    long double r_active = (long double)active_zone_radii[mean_idx];

                    mr1_active += r_active;
                    mr2_active += r_active*r_active;
                    mr3_active += r_active*r_active*r_active;
                    mr4_active += r_active*r_active*r_active*r_active;

                }

                mr1_active /= (long double)active_zone_radii_count;
                mr2_active /= (long double)active_zonerdii_count;
                mr3_active /= (long double)active_zone_radii_count;
                mr4_active /= (long double)active_zone_radii_count;

                M1_active[layer_idx] += mr1_active;
                M2_active[layer_idx] += mr2_active;
                M3_active[layer_idx] += mr3_active;
                M4_active[layer_idx] += mr4_active;

                varA_active = mr2_active - (mr1_active*mr1_active);
                k3A_active = mr3_active - 3*mr1_active*mr2_active + 2*(mr1_active*mr1_active*mr1_active);
                k4A_active = mr4_active - 4*mr1_active*mr3_active - 3*(mr2_active*mr2_active) + 12*(mr1_active*mr1_active)*mr2_active 
                                - 6*(mr1_active*mr1_active*mr1_active*mr1_active);

                varB_vector_active[layer_idx] += mr2_active - (mr1_active*mr1_active);
                k3B_vector_active[layer_idx] += mr3_active - 3*mr1_active*mr2_active + 2*(mr1_active*mr1_active*mr1_active);
                k4B_vector_active[layer_idx] += mr4_active - 4*mr1_active*mr3_active - 3*(mr2_active*mr2_active) 
                                            + 12*(mr1_active*mr1_active)*mr2_active - 6*(mr1_active*mr1_active*mr1_active*mr1_active);

                varA_vector_active[layer_idx] += varA_active;
                skewnessA_active[layer_idx] += k3A_active/(varA_active*sqrt(varA_active));
                kurtosisA_active[layer_idx] += k4A_active/(varA_active*varA_active);

                // Distribución de radio máximo

                rmax1 = max_r;
                rmax2 = max_r*max_r;
                rmax3 = max_r*max_r*max_r;
                rmax4 = max_r*max_r*max_r*max_r;

                M1_rmax[layer_idx] += rmax1;
                M2_rmax[layer_idx] += rmax2;
                M3_rmax[layer_idx] += rmax3;
                M4_rmax[layer_idx] += rmax4;

                varA_rmax = rmax2 - (rmax1*rmax1);
                k3A_rmax = rmax3 - 3*rmax1*rmax2 + 2*(rmax1*rmax1*rmax1);
                k4A_rmax = rmax4 - 4*rmax1*rmax3 - 3*(rmax2*rmax2) + 12*(rmax1*rmax1)*rmax2 - 6*(rmax1*rmax1*rmax1*rmax1);

                varB_vector_rmax[layer_idx] += rmax2 - (rmax1*rmax2);
                k3B_vector_rmax[layer_idx] += rmax3 - 3*rmax1*rmax2 + 2*(rmax1*rmax1*rmax1);
                k4B_vector_rmax[layer_idx] += rmax4 - 4*rmax1*rmax3 - 3*(rmax2*rmax2) + 12*(rmax1*rmax1)*rmax2 - 6*(rmax1*rmax1*rmax1*rmax1);

                varA_vector_rmax[layer_idx] += varA_rmax;
                skewnessA_rmax[layer_idx] += k3A_rmax/(varA_rmax*sqrt(varA_rmax));
                kurtosisA_rmax[layer_idx] += k4A_rmax/(varA_rmax*varA_rmax);

                // Todas las partículas

                long double ri_all = 0;
                long double mr1_all = 0, mr2_all = 0, mr3_all = 0, mr4_all = 0;

                // Iterar sobre TODAS las partículas del agregado

                for (int p_idx = 0; p_idx < particle_count; p_idx++) {
                    
                    // Simplemente lee el radio pre-calculado

                    ri_all = (long double)radii_list[p_idx];
                    
                    // Acumular las sumas para los momentos

                    mr1_all += ri_all;
                    mr2_all += ri_all*ri_all;
                    mr3_all += ri_all*ri_all*ri_all;
                    mr4_all += ri_all*ri_all*ri_all*ri_all;

                }

                mr1_all /= (long double)particle_count;
                mr2_all /= (long double)particle_count;
                mr3_all /= (long double)particle_count;
                mr4_all /= (long double)particle_count;

                M1_all[layer_idx] += mr1_all;
                M2_all[layer_idx] += mr2_all;
                M3_all[layer_idx] += mr3_all;
                M4_all[layer_idx] += mr4_all;

                long double varA_all = mr2_all - (mr1_all*mr1_all);
                long double k3A_all = mr3_all - 3*mr1_all*mr2_all + 2*(mr1_all*mr1_all*mr1_all);
                long double k4A_all = mr4_all - 4*mr1_all*mr3_all - 3*(mr2_all*mr2_all) 
                                    + 12*(mr1_all*mr1_all)*mr2_all - 6*(mr1_all*mr1_all*mr1_all*mr1_all);

                varA_vector_all[layer_idx] += varA_all;
                skewnessA_all[layer_idx] += k3A_all/(varA_all * sqrt(varA_all));
                kurtosisA_all[layer_idx] += k4A_all/(varA_all * varA_all);

                varB_vector_all[layer_idx] += mr2_all - (mr1_all*mr1_all);
                k3B_vector_all[layer_idx] += mr3_all - 3*mr1_all * mr2_all + 2*(mr1_all*mr1_all*mr1_all);
                k4B_vector_all[layer_idx] += mr4_all - 4*mr1_all * mr3_all - 3*(mr2_all*mr2_all) 
                                    + 12*(mr1_all*mr1_all)*mr2_all - 6*(mr1_all*mr1_all*mr1_all*mr1_all);

                // Radio de giro

                double current_rg = get_rg_sq(particle_count);
                total_rg_per_layer[layer_idx] += current_rg;

                // Acumular la media, varianza y el radio de giro

                // Calcula la métrica de esta capa para esta simulación m

                long double current_mean = mr1;
                long double current_rough = mr2 - (mr1*mr1);
                long double current_rmax = rmax2 - (rmax1*rmax1);
                double current_giro = get_rg_sq(particle_count);

                // Mean

                M1_mean_dist[layer_idx] += current_mean;
                M2_mean_dist[layer_idx] += current_mean*current_mean;
                M3_mean_dist[layer_idx] += current_mean*current_mean*current_mean;
                M4_mean_dist[layer_idx] += current_mean*current_mean*current_mean*current_mean;

                // Roughness

                M1_rough_dist[layer_idx] += current_rough;
                M2_rough_dist[layer_idx] += current_rough*current_rough;
                M3_rough_dist[layer_idx] += current_rough*current_rough*current_rough;
                M4_rough_dist[layer_idx] += current_rough*current_rough*current_rough*current_rough;

                // Giration radious

                M1_giro_dist[layer_idx] += current_giro;
                M2_giro_dist[layer_idx] += current_giro*current_giro;
                M3_giro_dist[layer_idx] += current_giro*current_giro*current_giro;
                M4_giro_dist[layer_idx] += current_giro*current_giro*current_giro*current_giro;

                // Max radius

                M1_rmax_dist[layer_idx] += current_rmax;
                M2_rmax_dist[layer_idx] += current_rmax*current_rmax;
                M3_rmax_dist[layer_idx] += current_rmax*current_rmax*current_rmax;
                M4_rmax_dist[layer_idx] += current_rmax*current_rmax*current_rmax*current_rmax;

                // Reiniciar el contador para la próxima capa de análisis

                active_zone_radii_count = 0;

            } // Fin if análisis

        } // Fin loop deposición

    } // Fin loop muestras


    char rg_filename[60];

    sprintf(rg_filename, "rg_dynamics_avg_T=%d_run%d.dat", (int)T, job_id);
    FILE* rg_file = fopen(rg_filename, "w");

    if (rg_file == NULL) {

        perror("No se pudo crear el archivo");
    } else {

        fprintf(rg_file, "# Particulas\tRg_promedio_RMS\n");

        for (int i = 0; i < layers; i++) {

            long long particles_at_layer = (i + 1) * PARTICLES_PER_LAYER;
            double avg_rg_sq = total_rg_per_layer[i] / N;
            double final_avg_rg = sqrt(avg_rg_sq);
            
            fprintf(rg_file, "%lld\t%f\n", particles_at_layer, final_avg_rg);
        }

        fclose(rg_file);
    }

    char file_roughnessC[60], file_skC[60], file_ktC[60], file_roughness_activeC[60],
    file_sk_activeC[60], file_kt_activeC[60], file_mean_r[60], file_mean_r_active[60], file_num_r[60],
    file_roughness_rmaxC[60], file_kt_rmaxC[60], file_sk_rmaxC[60], moments_interface[60], moments_active[60],
    moments_rmax[60], file_roughnessA[60], file_skA[60], file_ktA[60], file_roughnessB[60], file_skB[60], file_ktB[60],
    file_roughness_activeA[60], file_roughness_activeB[60], file_sk_activeA[60], file_sk_activeB[60], file_kt_activeA[60],
    file_kt_activeB[60], file_roughness_rmaxA[60], file_roughness_rmaxB[60], file_sk_rmaxA[60], file_sk_rmaxB[60],
    file_kt_rmaxA[60], file_kt_rmaxB[60], file_roughness_allA[60], file_roughness_allB[60], file_roughness_allC[60],
    file_sk_allA[60], file_sk_allB[60], file_sk_allC[60], file_kt_allA[60], file_kt_allB[60], file_kt_allC[60];

    char file_mean_mean[60], file_var_mean[60], file_mean_rough[60], file_var_rough[60], file_sk_mean[60], file_kt_mean[60],
    file_sk_rough[60], file_kt_rough[60], file_mean_giro[60], file_var_giro[60], file_sk_giro[60], file_kt_giro[60],
    file_mean_rmax[60], file_var_rmax[60], file_sk_rmax[60], file_kt_rmax[60];

    // Global

    sprintf(file_mean_r, "mean_r_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_mean_r_active, "mean_r_active_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_num_r, "num_per_layer_T=%d_run%d.dat", (int)T, job_id);
    sprintf(moments_interface, "moments_interfce_T=%d_run%d.dat", (int)T, job_id);
    sprintf(moments_active, "moments_active_T=%d_run%d.dat", (int)T, job_id);
    sprintf(moments_rmax, "moments_rmax_T=%d_run%d.dat", (int)T, job_id);

    // Método A

    sprintf(file_roughnessA, "roughness_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_activeA, "roughness_active_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_rmaxA, "roughness_rmax_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_allA, "roughness_all_A_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_skA, "skewness_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_activeA, "skewness_active_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_rmaxA, "skewness_rmax_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_allA, "skewness_all_A_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_ktA, "kurtosis_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_activeA, "kurtosis_active_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_rmaxA, "kurtosis_rmax_A_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_allA, "kurtosis_all_A_T=%d_run%d.dat", (int)T, job_id);

    // Método B

    sprintf(file_roughnessB, "roughness_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_activeB, "roughness_active_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_rmaxB, "roughness_rmax_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_allB, "roughness_all_B_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_skB, "skewness_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_activeB, "skewness_active_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_rmaxB, "skewness_rmax_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_allB, "skewness_all_B_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_ktB, "kurtosis_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_activeB, "kurtosis_active_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_rmaxB, "kurtosis_rmax_B_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_allB, "kurtosis_all_B_T=%d_run%d.dat", (int)T, job_id);

    // Método C

    sprintf(file_roughnessC, "roughness_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_activeC, "roughness_active_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_rmaxC, "roughness_rmax_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_roughness_allC, "roughness_all_C_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_skC, "skewness_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_activeC, "skewness_active_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_rmaxC, "skewness_rmax_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_allC, "skewness_all_C_T=%d_run%d.dat", (int)T, job_id);

    sprintf(file_ktC, "kurtosis_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_activeC, "kurtosis_active_C_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_rmaxC, "kurtosis_rmax_C_T=%d_run%d.dat", (int)T, job_id);    
    sprintf(file_kt_allC, "kurtosis_all_C_T=%d_run%d.dat", (int)T, job_id);

    // Distribución de la distribución (?)

    sprintf(file_mean_mean, "mean_mean_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_var_mean, "var_mean_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_mean, "skewness_mean_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_mean, "kurtosis_mean_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_mean_rough, "mean_roughness_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_var_rough, "var_roughness_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_rough, "skewness_rough_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_rough, "kurtosis_rough_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_mean_giro, "mean_giration_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_var_giro, "var_giration_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_giro, "skewness_giration_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_giro, "kurtosis_giration_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_mean_rmax, "mean_rmax_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_var_rmax, "var_rmax_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_sk_rmax, "sk_rmax_T=%d_run%d.dat", (int)T, job_id);
    sprintf(file_kt_rmax, "kt_rmax_T=%d_run%d.dat", (int)T, job_id);

    // Guardar archivos de método A y B

    FILE *filevarA = fopen(file_roughnessA, "w"), *fileskA = fopen(file_skA, "w"), *filektA = fopen(file_ktA, "w"), 
    *filevar_activeA = fopen(file_roughness_activeA, "w"), *filesk_activeA = fopen(file_sk_activeA, "w"),
    *filekt_activeA = fopen(file_kt_activeA, "w"), *filerrmaxA = fopen(file_roughness_rmaxA, "w"), 
    *fileskrmaxA = fopen(file_sk_rmaxA, "w"), *filektrmaxA = fopen(file_kt_rmaxA, "w"),
    *filevarB = fopen(file_roughnessB, "w"), *fileskB = fopen(file_skB, "w"), *filektB = fopen(file_ktB, "w"), 
    *filevar_activeB = fopen(file_roughness_activeB, "w"), *filesk_activeB = fopen(file_sk_activeB, "w"),
    *filekt_activeB = fopen(file_kt_activeB, "w"), *filerrmaxB = fopen(file_roughness_rmaxB, "w"), 
    *fileskrmaxB = fopen(file_sk_rmaxB, "w"), *filektrmaxB = fopen(file_kt_rmaxB, "w"), *filevarallA = fopen(file_roughness_allA, "w"),
    *filevarallB = fopen(file_roughness_allB, "w"), *fileskallA = fopen(file_sk_allA, "w"), *fileskallB = fopen(file_sk_allB, "w"),
    *filektallA = fopen(file_kt_allA, "w"), *filektallB = fopen(file_kt_allB, "w");

    for (int file_idx = 0; file_idx < layers; file_idx++){

        PARTICLES[file_idx] = (file_idx + 1)*PARTICLES_PER_LAYER;

        // Interfaz

        fprintf(fileskA, "%d %.9Lf\n", PARTICLES[file_idx], skewnessA[file_idx]/(long double)N);
        fprintf(filektA, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisA[file_idx]/(long double)N);
        fprintf(filevarA, "%d %.9Lf\n", PARTICLES[file_idx], varA_vector[file_idx]/(long double)N);

        varB = varB_vector[file_idx]/(long double)N;
        skewnessB = (k3B_vector[file_idx]/(long double)N)/(varB*sqrt(varB));
        kurtosisB = (k4B_vector[file_idx]/(long double)N)/(varB*varB);

        fprintf(filevarB, "%d %.9Lf\n", PARTICLES[file_idx], varB);
        fprintf(fileskB, "%d %.9Lf\n", PARTICLES[file_idx], skewnessB);
        fprintf(filektB, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisB);

        // Zona activa

        fprintf(filesk_activeA, "%d %.9Lf\n", PARTICLES[file_idx], skewnessA_active[file_idx]/(long double)N);
        fprintf(filekt_activeA, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisA_active[file_idx]/(long double)N);
        fprintf(filevar_activeA, "%d %.9Lf\n", PARTICLES[file_idx], varA_vector_active[file_idx]/(long double)N);

        varB_active = varB_vector_active[file_idx]/(long double)N;
        skewness_activeB = (k3B_vector_active[file_idx]/(long double)N)/(varB_active*sqrt(varB_active));
        kurtosis_activeB = (k4B_vector_active[file_idx]/(long double)N)/(varB_active*varB_active);

        fprintf(filevar_activeB, "%d %.9Lf\n", PARTICLES[file_idx], varB_active);
        fprintf(filesk_activeB, "%d %.9Lf\n", PARTICLES[file_idx], skewness_activeB);
        fprintf(filekt_activeB, "%d %.9Lf\n", PARTICLES[file_idx], kurtosis_activeB);

        // Radio maximo

        fprintf(fileskrmaxA, "%d %.9Lf\n", PARTICLES[file_idx], skewnessA_rmax[file_idx]/(long double)N);
        fprintf(filektrmaxA, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisA_rmax[file_idx]/(long double)N);
        fprintf(filerrmaxA, "%d %.9Lf\n", PARTICLES[file_idx], varA_vector_rmax[file_idx]/(long double)N);

        varB_rmax = varB_vector_rmax[file_idx]/(long double)N;
        skewness_rmaxB = (k3B_vector_rmax[file_idx]/(long double)N)/(varB_rmax*sqrt(varB_rmax));
        kurtosis_rmaxB = (k4B_vector_rmax[file_idx]/(long double)N)/(varB_rmax*varB_rmax);

        fprintf(filerrmaxB, "%d %.9Lf\n", PARTICLES[file_idx], varB_rmax);
        fprintf(fileskrmaxB, "%d %.9Lf\n", PARTICLES[file_idx], skewness_rmaxB);
        fprintf(filektrmaxB, "%d %.9Lf\n", PARTICLES[file_idx], kurtosis_rmaxB);

        // Todas las partículas

        fprintf(filevarallA, "%d %.9Lf\n", PARTICLES[file_idx], varA_vector_all[file_idx]/(long double)N);
        fprintf(fileskallA, "%d %.9Lf\n", PARTICLES[file_idx], skewnessA_all[file_idx]/(long double)N);
        fprintf(filektallA, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisA_all[file_idx]/(long double)N);

        long double varB_all = varB_vector_all[file_idx];
        long double skewness_allB = (k3B_vector_all[file_idx]/(long double)N)/(varB_all*sqrt(varB_all));
        long double kurtosis_allB = (k4B_vector_all[file_idx]/(long double)N)/(varB_all*varB_all);

        fprintf(filevarallB, "%d %.9Lf\n", PARTICLES[file_idx], varB_all);
        fprintf(fileskallB, "%d %.9Lf\n", PARTICLES[file_idx], skewness_allB);
        fprintf(filektallB, "%d %.9Lf\n", PARTICLES[file_idx], kurtosis_allB);

    }

    fclose(filevarA);
    fclose(fileskA);
    fclose(filektA);
    fclose(filevar_activeA);
    fclose(filesk_activeA);
    fclose(filekt_activeA);
    fclose(filerrmaxA);
    fclose(fileskrmaxA);
    fclose(filektrmaxA);
    fclose(filevarB);
    fclose(fileskB);
    fclose(filektB);
    fclose(filevar_activeB);
    fclose(filesk_activeB);
    fclose(filekt_activeB);
    fclose(fileskrmaxB);
    fclose(filektrmaxB);
    fclose(filevarallA);
    fclose(filevarallB);
    fclose(fileskallA);
    fclose(fileskallB);
    fclose(filektallA);
    fclose(filektallB);


    FILE *filevarC = fopen(file_roughnessC, "w"), *fileskC = fopen(file_skC, "w"), *filektC = fopen(file_ktC, "w"),
    *filevar_activeC = fopen(file_roughness_activeC, "w"), *filesk_activeC = fopen(file_sk_activeC, "w"),
    *filekt_activeC = fopen(file_kt_activeC, "w"), *filemean = fopen(file_mean_r, "w"), *filemean_act = fopen(file_mean_r_active, "w"),
    *filenumr = fopen(file_num_r, "w"), *filerrmaxC = fopen(file_roughness_rmaxC, "w"), *fileskrmaxC = fopen(file_sk_rmaxC, "w"),
    *filektrmaxC = fopen(file_kt_rmaxC, "w"), *filemint = fopen(moments_interface, "w"), *filemact = fopen(moments_active, "w"),
    *filemrmax = fopen(moments_rmax, "w"), *filevarallC = fopen(file_roughness_allC, "w"), *fileskallC = fopen(file_sk_allC, "w"),
    *filektallC = fopen(file_kt_allC, "w");

    for (int file_idx = 0; file_idx < layers; file_idx++){

        PARTICLES[file_idx] = (file_idx + 1)*PARTICLES_PER_LAYER;

        // Interface normalization

        M1[file_idx] /= (long double)N;
        M2[file_idx] /= (long double)N;
        M3[file_idx] /= (long double)N;
        M4[file_idx] /= (long double)N;

        // Rmax normalization

        M1_rmax[file_idx] /= (long double)N;
        M2_rmax[file_idx] /= (long double)N;
        M3_rmax[file_idx] /= (long double)N;
        M4_rmax[file_idx] /= (long double)N;

        // Active zone normalization

        M1_active[file_idx] /= (long double)N;
        M2_active[file_idx] /= (long double)N;
        M3_active[file_idx] /= (long double)N;
        M4_active[file_idx] /= (long double)N;

        // All particles normalization

        M1_all[file_idx] /= (long double)N;
        M2_all[file_idx] /= (long double)N;
        M3_all[file_idx] /= (long double)N;
        M4_all[file_idx] /= (long double)N;

        // Calcular la rugosidad con método C

        varC = M2[file_idx] - (M1[file_idx]*M1[file_idx]);

        k3C = M3[file_idx] - 3 * M1[file_idx] * M2[file_idx] + 2 * (M1[file_idx] * M1[file_idx] * M1[file_idx]);
        k4C = M4[file_idx] - 4 * M1[file_idx] * M3[file_idx] - 3 * M2[file_idx] * M2[file_idx] + 
                             12 * (M1[file_idx] * M1[file_idx]) * M2[file_idx] - 6 * (M1[file_idx] * M1[file_idx] * M1[file_idx] * M1[file_idx]);


        skewnessC = k3C / (varC * sqrt(varC));
        kurtosisC = k4C / (varC * varC);
        
        fprintf(filemean, "%d %.6Lf\n", PARTICLES[file_idx], M1[file_idx]);
        fprintf(fileskC, "%d %.9Lf\n", PARTICLES[file_idx], skewnessC);
        fprintf(filektC, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisC);
        fprintf(filevarC, "%d %.9Lf\n", PARTICLES[file_idx], varC);

        num_rt_particles[file_idx] /= (long double)N;

        // Calcular las cantidades estadísticas de la zona activa

        long double varC_active = M2_active[file_idx] - (M1_active[file_idx] * M1_active[file_idx]);
        
        long double k3C_active = M3_active[file_idx] - 3 * M1_active[file_idx] * M2_active[file_idx] + 2 * (M1_active[file_idx]*M1_active[file_idx]*M1_active[file_idx]);

        long double k4C_active = M4_active[file_idx] - 4 * M1_active[file_idx] * M3_active[file_idx] - 3 * (M2_active[file_idx]*M2_active[file_idx]) + 
                                                  12 * (M1_active[file_idx]*M1_active[file_idx]) * M2_active[file_idx] - 
                                                  6 * (M1_active[file_idx]*M1_active[file_idx]*M1_active[file_idx]*M1_active[file_idx]);

        long double skewnessC_active = k3C_active / (varC_active * sqrt(varC_active));
        long double kurtosisC_active = k4C_active / (varC_active * varC_active);

        fprintf(filemean_act, "%d %.6Lf\n", PARTICLES[file_idx], M1_active[file_idx]);
        fprintf(filevar_activeC, "%d %.9Lf\n", PARTICLES[file_idx], varC_active);
        fprintf(filesk_activeC, "%d %.9Lf\n", PARTICLES[file_idx], skewnessC_active);
        fprintf(filekt_activeC, "%d %.9Lf\n", PARTICLES[file_idx], kurtosisC_active);

        fprintf(filenumr, "%d %.6Lf\n", PARTICLES[file_idx], num_rt_particles[file_idx]/(long double)PARTICLES[file_idx]);

        // Cálculos para la distribución de radio máximo

        long double var_rmax = M2_rmax[file_idx] - (M1_rmax[file_idx]*M1_rmax[file_idx]);

        long double k3_rmax = M3_rmax[file_idx] - 3 * M1_rmax[file_idx] * M2_rmax[file_idx] + 
                                2*(M1_rmax[file_idx] * M1_rmax[file_idx] * M1_rmax[file_idx]);

        long double k4_rmax = M4_rmax[file_idx] - 4 * M1_rmax[file_idx] * M3_rmax[file_idx] - 3 * M2_rmax[file_idx] * M2_rmax[file_idx] + 
                                 12 * (M1_rmax[file_idx] * M1_rmax[file_idx]) * M2_rmax[file_idx] - 
                                 6 * (M1_rmax[file_idx] * M1_rmax[file_idx] * M1_rmax[file_idx] * M1_rmax[file_idx]);


        long double skewness_rmax = k3_rmax / (var_rmax * sqrt(var_rmax));
        long double kurtosis_rmax = k4_rmax / (var_rmax * var_rmax);

        fprintf(filerrmaxC, "%d %.6Lf\n", PARTICLES[file_idx], var_rmax);
        fprintf(fileskrmaxC, "%d %.9Lf\n", PARTICLES[file_idx], skewness_rmax);
        fprintf(filektrmaxC, "%d %.9Lf\n", PARTICLES[file_idx], kurtosis_rmax);

        // Cálculo para todas las partículas

        long double var_all = M2_all[file_idx] - (M1_all[file_idx]*M1_all[file_idx]);

        long double k3_all = M3_all[file_idx] - 3 * M1_all[file_idx] * M2_all[file_idx] + 
                                2*(M1_all[file_idx] * M1_all[file_idx] * M1_all[file_idx]);

        long double k4_all = M4_all[file_idx] - 4 * M1_all[file_idx] * M3_all[file_idx] - 3 * M2_all[file_idx] * M2_all[file_idx] + 
                                 12 * (M1_all[file_idx] * M1_all[file_idx]) * M2_all[file_idx] - 
                                 6 * (M1_all[file_idx] * M1_all[file_idx] * M1_all[file_idx] * M1_all[file_idx]);

        long double skewness_all = k3_all / (var_all * sqrt(var_all));
        long double kurtosis_all = k4_all / (var_all * var_all);

        fprintf(filevarallC, "%d %.6Lf\n", PARTICLES[file_idx], var_all);
        fprintf(fileskallC, "%d %.9Lf\n", PARTICLES[file_idx], skewness_all);
        fprintf(filektallC, "%d %.9Lf\n", PARTICLES[file_idx], kurtosis_all);

        // Guardar todos los momentos

        fprintf(filemint, "%d %.6Lf %.6Lf %.6Lf %.6Lf\n", PARTICLES[file_idx], M1[file_idx], M2[file_idx], M3[file_idx], M4[file_idx]);
        fprintf(filemact, "%d %.6Lf %.6Lf %.6Lf %.6Lf\n", PARTICLES[file_idx], M1_active[file_idx], M2_active[file_idx], M3_active[file_idx], M4_active[file_idx]);
        fprintf(filemrmax, "%d %.6Lf %.6Lf %.6Lf %.6Lf\n", PARTICLES[file_idx], M1_rmax[file_idx], M2_rmax[file_idx], M3_rmax[file_idx], M4_rmax[file_idx]);

    }

    fclose(filevarC);
    fclose(fileskC);
    fclose(filektC);
    fclose(filevar_activeC);
    fclose(filesk_activeC);
    fclose(filekt_activeC);
    fclose(filemean);
    fclose(filemean_act);
    fclose(filenumr);
    fclose(filerrmaxC);
    fclose(fileskrmaxC);
    fclose(filektrmaxC);
    fclose(filemint);
    fclose(filemact);
    fclose(filemrmax);
    fclose(filevarallC);
    fclose(fileskallC);
    fclose(filektallC);

    FILE *filemeanmean = fopen(file_mean_mean, "w"), *filevarmean = fopen(file_var_mean, "w"), *fileskmean = fopen(file_sk_mean, "w"),
    *filektmean = fopen(file_kt_mean, "w"), *filemeanrough = fopen(file_mean_rough, "w"), *filevarrough = fopen(file_var_rough, "w"),
    *fileskrough = fopen(file_sk_rough, "w"), *filektrough = fopen(file_kt_rough, "w"), *filemeangiro = fopen(file_mean_giro, "w"),
    *filevargiro = fopen(file_var_giro, "w"), *fileskgiro = fopen(file_sk_giro, "w"), *filektgiro = fopen(file_kt_giro, "w"),
    *filemeanrmax = fopen(file_mean_rmax, "w"), *filevarrmax = fopen(file_var_rmax, "w"), *fileskrmax = fopen(file_sk_rmax, "w"),
    *filektrmax = fopen(file_kt_rmax, "w");

    // Normalición de los datos

    for (int f = 0; f < layers; f++){

        // Número de partículas

        PARTICLES[f] = (f + 1)*PARTICLES_PER_LAYER;

        // Normalización de los datos

        M1_mean_dist[f] /= (long double)N;
        M2_mean_dist[f] /= (long double)N;
        M3_mean_dist[f] /= (long double)N;
        M4_mean_dist[f] /= (long double)N;

        M1_rough_dist[f] /= (long double)N;
        M2_rough_dist[f] /= (long double)N;
        M3_rough_dist[f] /= (long double)N;
        M4_rough_dist[f] /= (long double)N;

        M1_giro_dist[f] /= (long double)N;
        M2_giro_dist[f] /= (long double)N;
        M3_giro_dist[f] /= (long double)N;
        M4_giro_dist[f] /= (long double)N;

        M1_rmax_dist[f] /= (long double)N;
        M2_rmax_dist[f] /= (long double)N;
        M3_rmax_dist[f] /= (long double)N;
        M4_rmax_dist[f] /= (long double)N;

        // Realizar los cálculos de la varianza, cumulantes, skewness y kurtosis

        // Mean

        long double mean_mean = M1_mean_dist[f];

        long double var_mean = M2_mean_dist[f] - (M1_mean_dist[f]*M1_mean_dist[f]);

        long double k3_mean = M3_mean_dist[f] - 3*M1_mean_dist[f]*M2_mean_dist[f] + 2*(M1_mean_dist[f]*M1_mean_dist[f]*M1_mean_dist[f]);

        long double k4_mean = M4_mean_dist[f] - 4 * M1_mean_dist[f] * M3_mean_dist[f] - 3 * M2_mean_dist[f] * M2_mean_dist[f] + 
                             12 * (M1_mean_dist[f] * M1_mean_dist[f]) * M2_mean_dist[f] - 
                             6 * (M1_mean_dist[f] * M1_mean_dist[f] * M1_mean_dist[f] * M1_mean_dist[f]);


        long double skewness_mean = k3_mean / (var_mean*sqrt(var_mean));
        long double kurtosis_mean = k4_mean / (var_mean*var_mean);

        fprintf(filemeanmean, "%d %.9Lf\n", PARTICLES[f], mean_mean);
        fprintf(filevarmean, "%d %.9Lf\n", PARTICLES[f], var_mean);
        fprintf(fileskmean, "%d %.9Lf\n", PARTICLES[f], skewness_mean);
        fprintf(filektmean, "%d %.9Lf\n", PARTICLES[f], kurtosis_mean);


        // Roughness

        long double mean_rough = M1_rough_dist[f];

        long double var_rough = M2_rough_dist[f] - (M1_rough_dist[f]*M1_rough_dist[f]);

        long double k3_rough = M3_rough_dist[f] - 3*M1_rough_dist[f]*M2_rough_dist[f] + 2*(M1_rough_dist[f]*M1_rough_dist[f]*M1_rough_dist[f]);

        long double k4_rough = M4_rough_dist[f] - 4 * M1_rough_dist[f] * M3_rough_dist[f] - 3 * M2_rough_dist[f] * M2_rough_dist[f] + 
                             12 * (M1_rough_dist[f] * M1_rough_dist[f]) * M2_rough_dist[f] - 
                             6 * (M1_rough_dist[f] * M1_rough_dist[f] * M1_rough_dist[f] * M1_rough_dist[f]);


        long double skewness_rough = k3_rough / (var_rough*sqrt(var_rough));
        long double kurtosis_rough = k4_rough / (var_rough*var_rough);

        fprintf(filemeanrough, "%d %.9Lf\n", PARTICLES[f], mean_rough);
        fprintf(filevarrough, "%d %.9Lf\n", PARTICLES[f], var_rough);
        fprintf(fileskrough, "%d %.9Lf\n", PARTICLES[f], skewness_rough);
        fprintf(filektrough, "%d %.9Lf\n", PARTICLES[f], kurtosis_rough);

        // Radio de giro

        long double mean_giro = M1_giro_dist[f];

        long double var_giro = M2_giro_dist[f] - (M1_giro_dist[f]*M1_giro_dist[f]);

        long double k3_giro = M3_giro_dist[f] - 3*M1_giro_dist[f]*M2_giro_dist[f] + 2*(M1_giro_dist[f]*M1_giro_dist[f]*M1_giro_dist[f]);

        long double k4_giro = M4_giro_dist[f] - 4 * M1_giro_dist[f] * M3_giro_dist[f] - 3 * M2_giro_dist[f] * M2_giro_dist[f] + 
                             12 * (M1_giro_dist[f] * M1_giro_dist[f]) * M2_giro_dist[f] - 
                             6 * (M1_giro_dist[f] * M1_giro_dist[f] * M1_giro_dist[f] * M1_giro_dist[f]);


        long double skewness_giro = k3_giro / (var_giro*sqrt(var_giro));
        long double kurtosis_giro = k4_giro / (var_giro*var_giro);

        fprintf(filemeangiro, "%d %.9Lf\n", PARTICLES[f], mean_giro);
        fprintf(filevargiro, "%d %.9Lf\n", PARTICLES[f], var_giro);
        fprintf(fileskgiro, "%d %.9Lf\n", PARTICLES[f], skewness_giro);
        fprintf(filektgiro, "%d %.9Lf\n", PARTICLES[f], kurtosis_giro);

        // Radio máximo

        long double mean_r_max = M1_rmax_dist[f];

        long double var_r_max = M2_rmax_dist[f] - (M1_rmax_dist[f]*M1_rmax_dist[f]);

        long double k3_r_max = M3_rmax_dist[f] - 3*M1_rmax_dist[f]*M2_rmax_dist[f] + 2*(M1_rmax_dist[f]*M1_rmax_dist[f]*M1_rmax_dist[f]);

        long double k4_r_max = M4_rmax_dist[f] - 4 * M1_rmax_dist[f] * M3_rmax_dist[f] - 3 * M2_rmax_dist[f] * M2_rmax_dist[f] + 
                             12 * (M1_rmax_dist[f] * M1_rmax_dist[f]) * M2_rmax_dist[f] - 
                             6 * (M1_rmax_dist[f] * M1_rmax_dist[f] * M1_rmax_dist[f] * M1_rmax_dist[f]);

        long double skewness_r_max = k3_r_max / (var_r_max*sqrt(var_r_max));
        long double kurtosis_r_max = k4_r_max / (var_r_max*var_r_max);

        fprintf(filemeanrmax, "%d %.9Lf\n", PARTICLES[f], mean_r_max);
        fprintf(filevarrmax, "%d %.9Lf\n", PARTICLES[f], var_r_max);
        fprintf(fileskrmax, "%d %.9Lf\n", PARTICLES[f], skewness_r_max);
        fprintf(filektrmax, "%d %.9Lf\n", PARTICLES[f], kurtosis_r_max);

    }

    fclose(filemeanmean);
    fclose(filevarmean);
    fclose(fileskmean);
    fclose(filektmean);
    fclose(filemeanrough);
    fclose(filevarrough);
    fclose(fileskrough);
    fclose(filektrough);
    fclose(filemeangiro);
    fclose(filevargiro);
    fclose(fileskgiro);
    fclose(filektgiro);
    fclose(filemeanrmax);
    fclose(filevarrmax);
    fclose(fileskrmax);
    fclose(filektrmax);

    save_particles("particulas_final.dat");

    printf("\nSimulación completada.\n");

    return 0;
}
