import pygame
import random
import math
pygame.init()


max_density = 0
P = 100
class particle:
    mass = 1000
    radius = 5
    smooth_radius = 50
    density = 0
    preassure = [0,0]

    def __init__(self, x, y, vx, vy):
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy

    def speed_to_color(self):
        speed = (self.vx ** 2 + self.vy ** 2) ** 0.5
        t = max(0.0, min(0.001, speed))/2

        r = int(255 * t)
        g = 0
        b = int(255 * (1.0 - t))

        return (r, g, b)
    def density_to_color(self):
        t = min(1, self.density/max_density)

        r = int(255 * t)
        g = 0
        b = int(255 * (1.0 - t))

        return (r, g, b)

    def smoothingKernel(self, dst):
        if dst > self.smooth_radius:
            return 0
        val = self.smooth_radius - dst
        volume = 3.14159 * (self.smooth_radius ** 4) / 6
        return val * val / volume

    def smoothingKernelGradient(self, dst):
        if dst > self.smooth_radius:
            return 0
        val = dst-self.smooth_radius
        scale = 12/(3.14159 * (self.smooth_radius ** 4))
        return val * scale

    def calculateDensity(self, particles):
        self.density = 0
        for p in particles:
            dst = ((self.x - p.x) ** 2 + (self.y - p.y) ** 2) ** 0.5
            self.density += p.smoothingKernel(dst)*self.mass
        return self.density

    def calculatePreassure(self, particles):
        for p in particles:
            if p == self:
                continue
            dst = ((self.x - p.x) ** 2 + (self.y - p.y) ** 2) ** 0.5
            if dst == 0:
                continue
            grad = p.smoothingKernelGradient(dst)*self.mass
            dx = grad * (self.x - p.x)/dst
            dy = grad * (self.y - p.y)/dst
            self.preassure[0] -= dx/self.density
            self.preassure[1] -= dy/self.density
            p.preassure[0] += dx/p.density
            p.preassure[1] += dy/p.density




    def calculateForce(self, particles):
        self.vx += P*self.preassure[0]/self.density
        self.vy += P*self.preassure[1]/self.density
        self.vy += 9.8*dt


    def update(self, dt):
        self.x += self.vx * dt
        self.y += self.vy * dt

        if self.x < 0 or self.x > 1000:
            self.vx = -self.vx
        if self.y < 0 or self.y > 1000:
            self.vy = -self.vy
        self.preassure = [0,0]


    def draw(self,screen):
        pygame.draw.circle(screen, self.density_to_color(), (int(self.x), int(self.y)), self.radius)




def create_particles(n):
    particles = []
    for i in range(n):
        x = random.randint(0, 1000)
        y = random.randint(0, 1000)
        vx = random.uniform(-1, 1)
        vy = random.uniform(-1, 1)
        particles.append(particle(x, y, vx, vy))
    return particles


dt = 1
n = 1000
particles = create_particles(n)


screen = pygame.display.set_mode((1000, 1000))
running = True


while running:
    for evt in pygame.event.get():
        if evt.type == pygame.QUIT:
            running = False

    screen.fill((255, 255, 255))

    for p in particles:
        p.calculateDensity(particles)

    for p in particles:
        p.calculatePreassure(particles)
        p.calculateForce(particles)
        if p.density > max_density:
            max_density = p.density

    for p in particles:
        p.update(dt)
        p.draw(screen)
        start = (int(p.x), int(p.y))
        scale = 1000
        fx = p.preassure[0]
        fy = p.preassure[1]
        end = (int(p.x + fx * scale), int(p.y + fy * scale))
        pygame.draw.line(screen, (0,0,0), start, end, 2)
        ang = math.atan2(fy, fx)
        for sign in (+1, -1):
            theta = ang + sign * math.radians(30)
            dx = math.cos(theta) * 6
            dy = math.sin(theta) * 6
            pygame.draw.line(screen, (0,0,0), end, (end[0] - dx, end[1] - dy), 2)
    pygame.display.flip()

pygame.quit()
