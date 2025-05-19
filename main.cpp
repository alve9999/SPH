#include <SFML/Graphics.hpp>
#include <omp.h>
#include <iostream>
#include <cmath>
#define RADIUS 3.f

#define WIDTH 640
#define HEIGHT 480

#define SMOOTHING_RADIUS 20.f
#define MASS 1000.f

#define VISCOSITY_COEFFICIENT 2.0f
#define PREASSURE_COEFFICIENT 100.0f
#define MAX_VELOCITY 20.0f
#define DAMPING_COEFFICIENT 0.95f

using namespace std;

struct particle
{
    float x, y;
    float vx, vy;
    float density;
    float pressure_x;
    float pressure_y;
    float visc_x;
    float visc_y;
    int hash;
};


float smoothingKernel(float dst) {
    if (dst > SMOOTHING_RADIUS) return 0.f;
    float val    = SMOOTHING_RADIUS - dst;
    float volume = 3.14159f * pow(SMOOTHING_RADIUS,4) / 6.f;
    return val*val / volume;
}

float smoothingKernelGradient(float dst) {
    if (dst > SMOOTHING_RADIUS) return 0.f;
    float val   = dst - SMOOTHING_RADIUS;
    float scale = 12.f / (3.14159f * pow(SMOOTHING_RADIUS, 4));
    return val * scale;
}

void computeDensity(particle* P, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        particle& pi = P[i];
        pi.density = 0.f;
        for (int j = 0; j < N; j++) {
            particle& pj = P[j];
            float dx = pi.x - pj.x;
            float dy = pi.y - pj.y;
            float r  = std::sqrt(dx*dx + dy*dy);
            pi.density += MASS * smoothingKernel(r);
        }
    }
}

void computePreassure(particle* P, int N) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        particle& pi = P[i];
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            particle& pj = P[j];
            float dx = pi.x - pj.x;
            float dy = pi.y - pj.y;
            float r2 = dx*dx + dy*dy;
            if(r2 > SMOOTHING_RADIUS*SMOOTHING_RADIUS) continue;
            if(r2 < 0.0001f) continue;
            
            float r = std::sqrt(r2);
            float w = smoothingKernelGradient(r);
            
            // Pressure forces
            float pressure_force = PREASSURE_COEFFICIENT * (pi.density + pj.density) / 2.0f;
            pi.pressure_x -= MASS * w * (dx/r) * pressure_force / (pi.density * pj.density);
            pi.pressure_y -= MASS * w * (dy/r) * pressure_force / (pi.density * pj.density);
            
            // Viscosity forces (damping based on relative velocity)
            float dvx = pj.vx - pi.vx;
            float dvy = pj.vy - pi.vy;
            pi.visc_x += VISCOSITY_COEFFICIENT * MASS * (dvx) * smoothingKernel(r) / pj.density;
            pi.visc_y += VISCOSITY_COEFFICIENT * MASS * (dvy) * smoothingKernel(r) / pj.density;
        }
    }
}


sf::Color getParticleColor(const particle& p)
{
    float speed = sqrt(p.vx * p.vx + p.vy * p.vy);
    float t = min(1.f, max(0.f, speed/6.0f));
    return sf::Color(255 * t, 0, 255 * (1 - t));
}

void limitVelocity(float& vx, float& vy) {
    float v2 = vx*vx + vy*vy;
    if (v2 > MAX_VELOCITY*MAX_VELOCITY) {
        float scale = MAX_VELOCITY / sqrt(v2);
        vx *= scale;
        vy *= scale;
    }
}

void drawParticle(sf::RenderWindow& window, const particle& p)
{
    sf::CircleShape circle(RADIUS);
    circle.setFillColor(getParticleColor(p));
    circle.setPosition(p.x, p.y);
    window.draw(circle);
}

void updateParticle(particle& p, float dt)
{
    // Update velocity based on pressure
    p.vx += p.pressure_x * dt * PREASSURE_COEFFICIENT/p.density;
    p.vy += p.pressure_y * dt * PREASSURE_COEFFICIENT/p.density;
    p.vx *= 0.995f;
    p.vy *= 0.995f;
    
    limitVelocity(p.vx, p.vy);

    p.vy += 9.81f * dt;


    if (p.x < 0 || p.x > WIDTH) p.vx = -p.vx * DAMPING_COEFFICIENT;
    if (p.y < 0 || p.y > HEIGHT) p.vy = -p.vy * DAMPING_COEFFICIENT;

    p.x += p.vx * dt;
    p.y += p.vy * dt;
    
    p.hash = findParticleQuadrant(p.x, p.y);

    p.pressure_x = 0.f;
    p.pressure_y = 0.f;
    p.visc_x = 0.f;
    p.visc_y = 0.f;
}

int compare_by_hash(const void* a, const void* b) {
    int ha = ((particle*)a)->hash;
    int hb = ((particle*)b)->hash;
    return (ha > hb) - (ha < hb);
}

int findParticleQuadrant(float x, float y)
{
    int ix = (int)x;
    int iy = (int)y;
    int qx = ix % SMOOTHING_RADIUS;
    int qy = iy % SMOOTHING_RADIUS;
    return (qx * 2381 + qy * 6661) % (WIDTH/SMOOTHING_RADIUS * HEIGTH/SMOOTHING_RADIUS);
}

void initParticles(particle* particles, int N)
{
    for (int i = 0; i < N; ++i)
    {
        particles[i].x = WIDTH/4 + (rand() % (WIDTH/2));
        particles[i].y = HEIGHT/4 + (rand() % (HEIGHT/2));
        particles[i].vx = 0.f;
        particles[i].vy = 0.f;
        particles[i].pressure_x = 0.f;
        particles[i].pressure_y = 0.f;
        particles[i].density = 0.f;
        particles[i].visc_x = 0.f;
        particles[i].visc_y = 0.f;
        particles[i].hash = findParticleQuadrant(particles[i].x, particles[i].y);
    }
}


}

int main()
{
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Fluid");

    int N = 1000;
    float dt = 0.16f;

    
    int* lookup = (int*)malloc(sizeof(int) * (WIDTH/SMOOTHING_RADIUS * HEIGHT/SMOOTHING_RADIUS));
    particle* particles = (particle*)malloc(sizeof(particle) * N);
    initParticles(particles, N);
    sf::Clock frameClock;
    sf::Clock fpsClock;
    int    frameCount = 0;

    while (window.isOpen())
    {
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            sf::Vector2i mousePos = sf::Mouse::getPosition(window);
            for (int i = 0; i < N; i++) {
                float dx = mousePos.x - particles[i].x;
                float dy = mousePos.y - particles[i].y;
                float dist = sqrt(dx*dx + dy*dy);
                if (dist < 100.0f) {
                    float force = 30.0f;
                    particles[i].vx += dx * force * dt / dist;
                    particles[i].vy += dy * force * dt / dist;
                }
            }
        }
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        qsort(particles, N, sizeof(particle), compare_by_hash);

        computeDensity(particles, N);
        computePreassure(particles, N);

        for (int i = 0; i < N; i++)
        {
            updateParticle(particles[i], dt);
        }

        window.clear(sf::Color::White);
        for (int i = 0; i < N; i++)
        {
            drawParticle(window, particles[i]);
        }
        window.display();

        frameCount++;
        if (fpsClock.getElapsedTime().asSeconds() >= 1.0f) {
            float fps = frameCount / fpsClock.restart().asSeconds();
            std::cout << "FPS: " << fps << "\n";
            frameCount = 0;
        }
    }

    return 0;
}
