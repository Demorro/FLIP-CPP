#include <SFML/Graphics.hpp>
#include "Fluid.h"

#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/Font.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/Graphics/Text.hpp>
#include <SFML/System/Vector2.hpp>
#include <cstddef>
#include <string>
#include <algorithm>
#include <cmath>

// ---------------------------------------------------------------------------
// Configuration
// ---------------------------------------------------------------------------
static constexpr std::size_t GRID_WIDTH    = 18;
static constexpr std::size_t GRID_HEIGHT    = 18;
static constexpr std::size_t GRID_SPACING = 20;
static constexpr float       FIXED_DT  = 1.0f / 60.0f;
static constexpr unsigned    WINDOW_W = static_cast<unsigned>(GRID_WIDTH * GRID_SPACING);
static constexpr unsigned    WINDOW_H = static_cast<unsigned>(GRID_HEIGHT * GRID_SPACING);


int main()
{
    sf::RenderWindow window(
        sf::VideoMode(WINDOW_W, WINDOW_H),
        "FluidSim",
        sf::Style::Titlebar | sf::Style::Close
    );
    //window.setFramerateLimit(60);

    sim::Fluid<GRID_WIDTH, GRID_HEIGHT, GRID_SPACING> fluid{};

    sf::Clock clock;
    float accumulator = 0.0f;


     sf::Font font;
    font.loadFromFile("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf");

    sf::Text fpsText;
    fpsText.setFont(font);
    fpsText.setCharacterSize(14);
    fpsText.setFillColor(sf::Color::White);
    fpsText.setPosition(4.f, 4.f);
    sf::Clock fpsClock;
    int frameCount = 0;
    

    while (window.isOpen()) {
        // -- Events ----------------------------------------------------------
        sf::Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::KeyPressed &&
                event.key.code == sf::Keyboard::Escape)
                window.close();
            //if (event.type == sf::Event::KeyPressed &&
                //event.key.code == sf::Keyboard::Space)
                //pressures = fluid.step(FIXED_DT);
        }



        // -- Gravity from mouse position on click ------------------------------
        static constexpr float GRAVITY_STRENGTH = 600.0f;
        static sim::Vector2 gravity{0.0f, GRAVITY_STRENGTH}; // default: down
        sf::Vector2i mousePos = sf::Mouse::getPosition(window);
        float mx = mousePos.x - static_cast<float>(WINDOW_W) / 2.0f;
        float my = mousePos.y - static_cast<float>(WINDOW_H) / 2.0f;
        float len = std::sqrt(mx * mx + my * my);
        if (sf::Mouse::isButtonPressed(sf::Mouse::Left)) {
            if (len > 1.0f) {
                gravity = {(mx / len) * GRAVITY_STRENGTH, (my / len) * GRAVITY_STRENGTH};
            }
        }

        // -- Update (fixed timestep) -----------------------------------------
        accumulator += clock.restart().asSeconds();
        while (accumulator >= FIXED_DT) {
            fluid.Step(FIXED_DT, gravity);
            accumulator -= FIXED_DT;
        }

        // -- Render ----------------------------------------------------------
        window.clear(sf::Color::Black);

        // Grid cells — one batched draw instead of 6400 individual shapes
        sf::VertexArray cellVerts(sf::Quads, GRID_WIDTH * GRID_HEIGHT * 4);
        const std::array<std::array<float, GRID_HEIGHT>, GRID_WIDTH>& pressures = fluid.GetPressures();
        const float cellSize = static_cast<float>(GRID_SPACING) - 1.0f;
        for (std::size_t x = 0; x < GRID_WIDTH; ++x) {
            for (std::size_t y = 0; y < GRID_HEIGHT; ++y) {
                const std::size_t i = (x * GRID_HEIGHT + y) * 4;
                const float px = static_cast<float>(x * GRID_SPACING);
                const float py = static_cast<float>(y * GRID_SPACING);

                cellVerts[i + 0].position = {px, py};
                cellVerts[i + 1].position = {px + cellSize, py};
                cellVerts[i + 2].position = {px + cellSize, py + cellSize};
                cellVerts[i + 3].position = {px, py + cellSize};

                sf::Color color;
                if (fluid.GetCell(x, y).GetCellState() == sim::CellState::Solid) {
                    color = sf::Color::White;
                } else if (fluid.GetCell(x, y).GetCellState() == sim::CellState::Fluid) {
                    float p = pressures[x][y];
                    static constexpr float maxP = 1500000.0f;
                    float t = std::min(std::max(p / maxP, -1.0f), 1.0f);
                    uint8_t r = static_cast<uint8_t>(std::max(t, 0.0f) * 255);
                    uint8_t b = static_cast<uint8_t>(std::max(-t, 0.0f) * 255);
                    //color = sf::Color(r, 0, b);
                    color = sf::Color::Red;
                }
                else {
                    color = sf::Color::Black;
                }
                cellVerts[i + 0].color = color;
                cellVerts[i + 1].color = color;
                cellVerts[i + 2].color = color;
                cellVerts[i + 3].color = color;
            }
        }
        window.draw(cellVerts);

        // Particles — batched as small quads
        const auto& particles = fluid.GetParticles();
        sf::VertexArray particleVerts(sf::Quads, particles.size() * 4);
        static constexpr float DOT_RADIUS = 8.0f;
        for (std::size_t i = 0; i < particles.size(); ++i) {
            const float px = particles[i].position.x;
            const float py = particles[i].position.y;
            const std::size_t vi = i * 4;
            particleVerts[vi + 0].position = {px - DOT_RADIUS, py - DOT_RADIUS};
            particleVerts[vi + 1].position = {px + DOT_RADIUS, py - DOT_RADIUS};
            particleVerts[vi + 2].position = {px + DOT_RADIUS, py + DOT_RADIUS};
            particleVerts[vi + 3].position = {px, py + DOT_RADIUS};
            particleVerts[vi + 0].color = sf::Color::Cyan;
            particleVerts[vi + 1].color = sf::Color::Cyan;
            particleVerts[vi + 2].color = sf::Color::Cyan;
            particleVerts[vi + 3].color = sf::Color::Cyan;
        }
        //window.draw(particleVerts);

        // -- Gravity direction arrow -----------------------------------------
        if (len > 1.0f) {
            const float cx = static_cast<float>(WINDOW_W) / 2.0f;
            const float cy = static_cast<float>(WINDOW_H) / 2.0f;
            const float arrowLen = 40.0f;
            const float gx = (gravity.x / GRAVITY_STRENGTH) * arrowLen;
            const float gy = (gravity.y / GRAVITY_STRENGTH) * arrowLen;

            // Shaft
            sf::VertexArray shaft(sf::Lines, 2);
            shaft[0].position = {cx, cy};
            shaft[1].position = {cx + gx, cy + gy};
            shaft[0].color = sf::Color::Yellow;
            shaft[1].color = sf::Color::Yellow;
            window.draw(shaft);

            // Arrowhead — two small lines angled off the tip
            const float headLen = 10.0f;
            const float nx = gx / arrowLen;
            const float ny = gy / arrowLen;
            // Rotate ±135 degrees from the direction vector
            const float cos135 = -0.7071f;
            const float sin135 = 0.7071f;
            float h1x = (nx * cos135 - ny * sin135) * headLen;
            float h1y = (nx * sin135 + ny * cos135) * headLen;
            float h2x = (nx * cos135 + ny * sin135) * headLen;
            float h2y = (-nx * sin135 + ny * cos135) * headLen;

            sf::VertexArray head(sf::Lines, 4);
            head[0].position = {cx + gx, cy + gy};
            head[1].position = {cx + gx + h1x, cy + gy + h1y};
            head[2].position = {cx + gx, cy + gy};
            head[3].position = {cx + gx + h2x, cy + gy + h2y};
            head[0].color = sf::Color::Yellow;
            head[1].color = sf::Color::Yellow;
            head[2].color = sf::Color::Yellow;
            head[3].color = sf::Color::Yellow;
            window.draw(head);
        }

         // -- FPS counter -----------------------------------------------------
        ++frameCount;
        if (fpsClock.getElapsedTime().asSeconds() >= 1.0f) {
            fpsText.setString("FPS: " + std::to_string(frameCount));
            frameCount = 0;
            fpsClock.restart();
        }
        window.draw(fpsText);

        window.display();
    }

    return 0;
}
