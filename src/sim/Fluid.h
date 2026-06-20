#pragma once

#include <cmath>
#include <cstddef>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <utility>
#include <algorithm>
#include <sstream>
#include <limits>
#include <vector>
#include <memory>


namespace sim {
struct Vector2{
    float x;
    float y;

    Vector2 Normalized() const { float len = std::sqrt(x * x + y * y); return {x / len, y / len}; }
};

inline Vector2 operator+(Vector2 a, Vector2 b) { return {a.x + b.x, a.y + b.y}; }
inline Vector2 operator-(Vector2 a, Vector2 b) { return {a.x - b.x, a.y - b.y}; }
inline Vector2 operator*(Vector2 a, float s)   { return {a.x * s, a.y * s}; }
inline Vector2 operator*(float s, Vector2 a)   { return {a.x * s, a.y * s}; }
inline Vector2 operator/(Vector2 a, float s)   { return {a.x / s, a.y / s}; }
inline float Length(Vector2 a)                  { return std::sqrt(a.x * a.x + a.y * a.y); }
inline float Distance(Vector2 a, Vector2 b)    { return Length(a - b); }
inline float DistanceSqr(Vector2 a, Vector2 b) { float dx = a.x - b.x; float dy = a.y - b.y; return dx * dx + dy * dy; }
}


/*
 * ////////////////////////////////////////////////////////////////////////////
 * ///// PARTICLE //////
 * ////////////////////////////////////////////////////////////////////////////
 */
namespace sim::particle
{
    struct Particle{
    Vector2 position{0,0};
    Vector2 velocity{0,0};

    Particle() = default;
    Particle(Vector2 pos, Vector2 vel) : position(pos), velocity(vel) {}

    void Simulate(float dt, Vector2 gravity, Vector2 domainMin, Vector2 domainMax){
        /* Semi-implicit Euler integration, if you want to be all fancy about it. */
        velocity = velocity + (dt * gravity);
        position = position + (dt * velocity);
        position.x = std::clamp(position.x, domainMin.x, domainMax.x);
        position.y = std::clamp(position.y, domainMin.y, domainMax.y);
    }

    void Move(Vector2 motion) {
        position = position + motion;
    }

};
}

/*
 * ////////////////////////////////////////////////////////////////////////////
 * ///// CELL //////
 * ////////////////////////////////////////////////////////////////////////////
 */
namespace sim::cell{
/*
 * Cell type to provide utility methods for accessing a single cell.
 * References mutable velocities buffers, does not store them, as velocities are shared between grid boundaries.
 */

enum class Direction 
{
    UP,
    DOWN,
    LEFT,
    RIGHT
};

inline Direction OppositeDirection(Direction dir){
    switch (dir) {
        case Direction::UP:
            return Direction::DOWN;
        case Direction::DOWN:
            return Direction::UP;
        case Direction::LEFT:
            return Direction::RIGHT;
        case Direction::RIGHT:
            return Direction::LEFT;
        default:
            throw std::runtime_error("Impossible case");
    }
}

inline Direction RotateDirectionClockwise(Direction dir){
    switch (dir) {
        case Direction::UP:
            return Direction::RIGHT;
        case Direction::DOWN:
            return Direction::LEFT;
        case Direction::LEFT:
            return Direction::UP;
        case Direction::RIGHT:
            return Direction::DOWN;
        default:
            throw std::runtime_error("Impossible case");
    }
}

enum class CellState
{
    Solid,
    Fluid,
    Air
};

struct CellVelocities{
    float& leftVelocity;
    float& rightVelocity;
    float& upVelocity;
    float& downVelocity;

    float Divergence() const {
        return rightVelocity - leftVelocity + downVelocity - upVelocity;
    }
};

/* The same data as CellVelocities, but this is used to get u1, u2, u3 and u4 of the staggered grid
   These ones are pointers, because it is possible for these values to be undefined 
   We also provide deltaX and deltaY, because we're working out deltas to figure these out anyway,
   and we need them to do particle->grid velocity transfer (to calculate the weights) */
class StaggeredCornerVelocities{ 
public:
    StaggeredCornerVelocities(float* corner1, size_t corner1FlatIndex, float* corner2, size_t corner2FlatIndex, float* corner3, size_t corner3FlatIndex, float* corner4, size_t corner4FlatIndex, float deltaX, float deltaY):
        corner1(corner1),
        corner2(corner2),
        corner3(corner3),
        corner4(corner4),
        corner1FlatIndex(corner1FlatIndex),
        corner2FlatIndex(corner2FlatIndex),
        corner3FlatIndex(corner3FlatIndex),
        corner4FlatIndex(corner4FlatIndex),
        deltaX(deltaX),
        deltaY(deltaY)
    {
    }

    float* corner1 = nullptr;
    float* corner2 = nullptr;
    float* corner3 = nullptr;
    float* corner4 = nullptr;
    size_t corner1FlatIndex = std::numeric_limits<size_t>::max();
    size_t corner2FlatIndex = std::numeric_limits<size_t>::max();
    size_t corner3FlatIndex = std::numeric_limits<size_t>::max();
    size_t corner4FlatIndex = std::numeric_limits<size_t>::max();

    float deltaX = 0;
    float deltaY = 0;

    uint8_t NumDefinedVels() const {
        uint8_t counter = 0;
        if (corner1 != nullptr){++counter;}
        if (corner2 != nullptr){++counter;}
        if (corner3 != nullptr){++counter;}
        if (corner4 != nullptr){++counter;}
        return counter;
    }
};

template <size_t GRID_CELL_SIZE>
class Cell {
public:

    Cell(size_t xIndex, size_t yIndex):
        m_xIndex(xIndex),
        m_yIndex(yIndex)
        {
            m_leftGlobalPos = m_xIndex * GRID_CELL_SIZE;
            m_topGlobalPos = m_yIndex * GRID_CELL_SIZE;
            m_rightGlobalPos = m_leftGlobalPos + GRID_CELL_SIZE;
            m_bottomGlobalPos = m_topGlobalPos + GRID_CELL_SIZE;
            m_centerPos = Vector2{m_leftGlobalPos + static_cast<float>(GRID_CELL_SIZE)/2, m_topGlobalPos + static_cast<float>(GRID_CELL_SIZE)/2};
        }

    void HookupAdjacents(Cell* leftCell, Cell* rightCell, Cell* upCell, Cell* downCell) {
        m_leftCell = leftCell;
        m_rightCell = rightCell;
        m_upCell = upCell;
        m_downCell = downCell;
    }

    bool IsInsideCell(Vector2 globalPoint) const{

        if (globalPoint.x >= m_leftGlobalPos && globalPoint.x < m_rightGlobalPos){
            if (globalPoint.y >= m_topGlobalPos && globalPoint.y < m_bottomGlobalPos){
                return true;
            }
        }

        return false;
    }

    size_t GetXIndex() const {return m_xIndex;}
    size_t GetYIndex() const {return m_yIndex;}

    bool IsLeftSurfaceCompressible() const { return m_leftCell != nullptr && m_leftCell->GetCellState() != CellState::Solid; }
    bool IsRightSurfaceCompressible() const { return m_rightCell != nullptr && m_rightCell->GetCellState() != CellState::Solid; }
    bool IsUpSurfaceCompressible() const { return m_upCell != nullptr && m_upCell->GetCellState() != CellState::Solid; }
    bool IsDownSurfaceCompressible() const { return m_downCell != nullptr && m_downCell->GetCellState() != CellState::Solid; }
    bool IsSurfaceCompresssible(Direction dir) const {
        switch (dir) {
            case Direction::LEFT:
                return IsLeftSurfaceCompressible();
            case Direction::RIGHT:
                return IsRightSurfaceCompressible();
            case Direction::UP:
                return IsUpSurfaceCompressible();
            case Direction::DOWN: 
                return IsDownSurfaceCompressible();
            default:
                throw std::runtime_error("Impossible case");
        }
    }

    size_t NumCompressibleSurfaces() const {
        size_t counter = 0;
        if (IsLeftSurfaceCompressible()){++counter;}
        if (IsRightSurfaceCompressible()){++counter;}
        if (IsUpSurfaceCompressible()){++counter;}
        if (IsDownSurfaceCompressible()){++counter;}
        return counter;
    }

    /*
     * Get a X/Y vector with the distances to the closest boundaries in this cell, a
     * along with a pair of directions telling us what the closes boundaries are.
     * The first element of the pair is Left/Right, the second Up/Down
     */
    std::pair<Vector2, std::pair<Direction, Direction>> DistancesFromClosestBoundaries(Vector2 globalPointInsideCell) const {
        if (!IsInsideCell(globalPointInsideCell)){
            std::ostringstream oss;
            oss << "globalPointInsideCell is not inside cell XPos: " << globalPointInsideCell.x << " YPos: " << globalPointInsideCell.y << " xIndex: " << m_xIndex << " yIndex: " << m_yIndex;
            throw std::runtime_error(oss.str());
        }

        Direction leftOrRight = globalPointInsideCell.x < m_centerPos.x ? Direction::LEFT : Direction::RIGHT;
        Direction upOrDown = globalPointInsideCell.y < m_centerPos.y ? Direction::UP : Direction::DOWN;

        const float deltaX = leftOrRight == Direction::LEFT ? std::abs(globalPointInsideCell.x - m_leftGlobalPos) :
                                                              std::abs(globalPointInsideCell.x - m_rightGlobalPos);
        
        const float deltaY = upOrDown == Direction::UP ? std::abs(globalPointInsideCell.y - m_topGlobalPos) :
                                                         std::abs(globalPointInsideCell.y - m_bottomGlobalPos);

        return std::make_pair(Vector2{deltaX, deltaY}, std::make_pair(leftOrRight, upOrDown));
    }

    const Cell& GetAdjacentCell(Direction direction) const {
        switch (direction) {
            case Direction::UP:
                return GetAdjacentCell(0, -1);
            case Direction::DOWN:
                return GetAdjacentCell(0, 1);
            case Direction::LEFT:
                return GetAdjacentCell(-1, 0);  
            case Direction::RIGHT:
                return GetAdjacentCell(1, 0);
            default:
                throw std::runtime_error("Impossible case");
        }
    }

    const Cell& GetAdjacentCell(int16_t xOffset, int16_t yOffset) const{
        const Cell* cellPtr = this;

        while (xOffset != 0){ 
            if (xOffset < 0){
                cellPtr = cellPtr->m_leftCell;
                if (cellPtr == nullptr){
                    throw std::runtime_error("Error accessing left cell");
                }
                xOffset += 1;
            }
            else if (xOffset > 0){
                cellPtr = cellPtr->m_rightCell;
                if (cellPtr == nullptr){
                    throw std::runtime_error("Error accessing right cell");
                }
                xOffset -= 1;
            }      
        }

        while (yOffset != 0){ 
            if (yOffset < 0){
                cellPtr = cellPtr->m_upCell;
                if (cellPtr == nullptr){
                    throw std::runtime_error("Error accessing up cell");
                }
                yOffset += 1;
            }
            else if (yOffset > 0){
                cellPtr = cellPtr->m_downCell;
                if (cellPtr == nullptr){
                    throw std::runtime_error("Error accessing down cell");
                }
                yOffset -= 1;
            }      
        }

        return *cellPtr;
    }


    void SetCellState(CellState state){ m_cellState = state;}
    CellState GetCellState() const { return m_cellState; }

    bool IsAirCell(const std::vector<sim::particle::Particle>& particles){
        /* Really bad performance wise, fix later */
        for(const auto& particle : particles){
            if (IsInsideCell(particle.position)){
                return false;
            }
        }
        return true;
    }

private:
    Cell* m_leftCell = nullptr;
    Cell* m_rightCell = nullptr;
    Cell* m_upCell = nullptr;
    Cell* m_downCell = nullptr;
    size_t m_xIndex = 0;
    size_t m_yIndex = 0;
    size_t m_leftGlobalPos = 0;
    size_t m_topGlobalPos = 0;
    size_t m_rightGlobalPos = 0;
    size_t m_bottomGlobalPos = 0;
    Vector2 m_centerPos{0,0};

    CellState m_cellState = CellState::Fluid;
    
};
}

namespace sim {

    using namespace sim::cell;
    using namespace sim::particle;
    

    template <size_t GRID_WIDTH, size_t GRID_HEIGHT, size_t GRID_CELL_SIZE>
    class Fluid
    {
    private:
        bool ValidCell(int16_t x, int16_t y){
            // size_t makes the >= 0 checks redundant, it just feels _wrong_ without them though I dunno.
            if (x >= 0 && x < m_cells.size()){
                if (y >= 0 && y < m_cells[x].size()){
                    return true;
                }
            }
            return false;
        }
    public:
        Fluid(){
            InitCells();
            HookupCellConnections();
            SetSolidBorders();
            SeedParticles();
        }

        void Step(float dt, Vector2 gravity){
            SimulateParticles(dt, gravity);
            ClassifyCells();
            PushApartParticles();
            VelocityTransferParticlesToGrid();
            EnforceBoundaryConditions();
            CalculateGridDensities();


            //Save old velocities for FLIP method
            std::array<float, (GRID_WIDTH + 1) * GRID_HEIGHT> oldHorizontalVelocities = m_horizontalVelocities;
            std::array<float, GRID_WIDTH * (GRID_HEIGHT + 1)> oldVerticalVelocities = m_verticalVelocities;

            MakeIncompressible(dt);

            VelocityTransferGridtoParticles(oldHorizontalVelocities, oldVerticalVelocities);
            PushApartParticles();
        }

        const Cell<GRID_CELL_SIZE>& GetCellFromGlobalPos(Vector2 globalPos){

            size_t xIndex = globalPos.x / GRID_CELL_SIZE;
            size_t yIndex = globalPos.y / GRID_CELL_SIZE;

            if(!ValidCell(xIndex, yIndex)){
                throw std::runtime_error("Attempting to access invalid cell");
            }

            return *m_cells[xIndex][yIndex].get();
        }


        const Cell<GRID_CELL_SIZE>& GetCell(int16_t xIndex, int16_t yIndex){
            if(!ValidCell(xIndex, yIndex)){
                throw std::runtime_error("Attempting to access invalid cell");
            }

            return *m_cells[xIndex][yIndex].get();
        }

        const std::vector<Particle>& GetParticles() const { return m_particles; }

        const std::array<std::array<float, GRID_HEIGHT>, GRID_WIDTH>& GetPressures() { return m_pressures;}
    private:

        void ClassifyCells() {
            // Reset all non-solid cells to Air
            for (size_t x = 0; x < GRID_WIDTH; ++x) {
                for (size_t y = 0; y < GRID_HEIGHT; ++y) {
                    Cell<GRID_CELL_SIZE>& cell = *m_cells[x][y];
                    if (cell.GetCellState() != CellState::Solid)
                        cell.SetCellState(CellState::Air);
                }
            }

            m_particlesInCells = {};

            // Mark cells containing particles as Fluid
            for (auto& particle : m_particles) {
                size_t cellXIndex = static_cast<size_t>(particle.position.x / GRID_CELL_SIZE);
                size_t cellYIndex = static_cast<size_t>(particle.position.y / GRID_CELL_SIZE);

                m_particlesInCells[cellXIndex][cellYIndex].push_back(&particle);

                if (ValidCell(cellXIndex, cellYIndex) && m_cells[cellXIndex][cellYIndex]->GetCellState() != CellState::Solid){
                    m_cells[cellXIndex][cellYIndex]->SetCellState(CellState::Fluid);
                }
            }
        }

        void PushTwoParticlesApart(Particle& particleA, Particle& particleB, const float PUSH_APART_DIST_SQ, const float PUSH_APART_DISTANCE){
            const float dx = particleA.position.x - particleB.position.x;
            const float dy = particleA.position.y - particleB.position.y;
            const float distSq = dx * dx + dy * dy;
            if (distSq < PUSH_APART_DIST_SQ && distSq > 0.0f){
                const float dist = std::sqrt(distSq);
                const float overlap = (PUSH_APART_DISTANCE - dist) / 3.0f;
                const float scale = overlap / dist;
                particleA.position.x += dx * scale;
                particleA.position.y += dy * scale;
                particleB.position.x -= dx * scale;
                particleB.position.y -= dy * scale;
            }
        }

        //Call after ClassifyCells such that m_particleSInCells is populated
        void PushApartParticles(){
            static const std::vector<Particle*> empty;
            for (size_t x = 0; x < GRID_WIDTH; ++x) {
                for (size_t y = 0; y < GRID_HEIGHT; ++y) {
                    const std::vector<Particle*>& thisCell       = m_particlesInCells[x][y];
                    const std::vector<Particle*>& rightCell      = x + 1 < GRID_WIDTH  ? m_particlesInCells[x+1][y]   : empty;
                    const std::vector<Particle*>& belowCell      = y + 1 < GRID_HEIGHT ? m_particlesInCells[x][y+1]   : empty;
                    const std::vector<Particle*>& belowLeftCell  = x > 0 && y + 1 < GRID_HEIGHT ? m_particlesInCells[x-1][y+1] : empty;
                    const std::vector<Particle*>& belowRightCell = x + 1 < GRID_WIDTH && y + 1 < GRID_HEIGHT ? m_particlesInCells[x+1][y+1] : empty;

                    static constexpr float PUSH_APART_DISTANCE = GRID_CELL_SIZE/2.5f;
                    static constexpr float PUSH_APART_DIST_SQ = PUSH_APART_DISTANCE * PUSH_APART_DISTANCE;
                    for (size_t i = 0; i < thisCell.size(); ++i) {
                        Particle* particle = thisCell[i];

                        //push against other particles in thisCell (only j > i to avoid duplicate pairs)
                        for (size_t j = i + 1; j < thisCell.size(); ++j) {
                            PushTwoParticlesApart(*particle, *thisCell[j], PUSH_APART_DIST_SQ, PUSH_APART_DISTANCE);
                        }
                        for (size_t j = 0; j < rightCell.size(); ++j) {
                            PushTwoParticlesApart(*particle, *rightCell[j], PUSH_APART_DIST_SQ, PUSH_APART_DISTANCE);
                        }
                        for (size_t j = 0; j < belowCell.size(); ++j) {
                            PushTwoParticlesApart(*particle, *belowCell[j], PUSH_APART_DIST_SQ, PUSH_APART_DISTANCE);
                        }
                        for (size_t j = 0; j < belowLeftCell.size(); ++j) {
                            PushTwoParticlesApart(*particle, *belowLeftCell[j], PUSH_APART_DIST_SQ, PUSH_APART_DISTANCE);
                        }
                        for (size_t j = 0; j < belowRightCell.size(); ++j) {
                            PushTwoParticlesApart(*particle, *belowRightCell[j], PUSH_APART_DIST_SQ, PUSH_APART_DISTANCE);
                        }

                    }
                }
            }

            // Clamp particles back into the fluid domain in case pushing moved them into walls
            static constexpr float minBound = static_cast<float>(GRID_CELL_SIZE);
            static constexpr float maxBoundX = static_cast<float>((GRID_WIDTH - 1) * GRID_CELL_SIZE);
            static constexpr float maxBoundY = static_cast<float>((GRID_HEIGHT - 1) * GRID_CELL_SIZE);
            for (auto& particle : m_particles) {
                particle.position.x = std::clamp(particle.position.x, minBound, maxBoundX);
                particle.position.y = std::clamp(particle.position.y, minBound, maxBoundY);
            }
        }

        void CalculateGridDensities(){
            /* All velocity based solvers are vulnerable to drift.
             * This attempts to correct that by applying outwards expansion forces until
             * the cell densities are similar to the densities when the simulation began  */

            m_densities = {};
            //Accumulate weighted velocities from particles
            for (const auto& particle : m_particles){
                const Cell<GRID_CELL_SIZE>& cell = GetCellFromGlobalPos(particle.position);

                // Work out which cell we're going to need to be reaching into, it's the closest
                std::pair<Vector2, std::pair<Direction, Direction>> directionData = cell.DistancesFromClosestBoundaries(particle.position);

                const Vector2& edgeDistances = directionData.first;   // (dx_to_edge, dy_to_edge), both in [0, h/2]
                const Direction horizontalDir = directionData.second.first;
                const Direction verticalDir   = directionData.second.second;

                const float h = static_cast<float>(GRID_CELL_SIZE);

                const float tx = (horizontalDir == Direction::RIGHT)
                    ? (0.5f - edgeDistances.x / h)   // particle in right half: tx ∈ [0, 0.5]
                    : (0.5f + edgeDistances.x / h);  // particle in left half:  tx ∈ [0.5, 1]

                const float ty = (verticalDir == Direction::DOWN)
                    ? (0.5f - edgeDistances.y / h)   // particle in lower half: ty ∈ [0, 0.5]
                    : (0.5f + edgeDistances.y / h);  // particle in upper half: ty ∈ [0.5, 1]

                const float sx = 1.0f - tx;
                const float sy = 1.0f - ty;

                const float w1 = sx * sy; 
                const float w2 = tx * sy; 
                const float w3 = tx * ty; 
                const float w4 = sx * ty; 

                const Cell<GRID_CELL_SIZE>& verticalCell = cell.GetAdjacentCell(verticalDir);
                const Cell<GRID_CELL_SIZE>& horizontalCell = cell.GetAdjacentCell(horizontalDir);

                const int x0Index = (horizontalDir == Direction::RIGHT) ? cell.GetXIndex() : horizontalCell.GetXIndex();
                const int x1Index = (horizontalDir == Direction::RIGHT) ? horizontalCell.GetXIndex() : cell.GetXIndex();

                const int y0Index = (verticalDir == Direction::DOWN)    ? cell.GetYIndex() : verticalCell.GetYIndex();
                const int y1Index = (verticalDir == Direction::DOWN)    ? verticalCell.GetYIndex() : cell.GetYIndex();

                m_densities[x0Index][y0Index] += w1;
                m_densities[x1Index][y0Index] += w2;
                m_densities[x1Index][y1Index] += w3;
                m_densities[x0Index][y1Index] += w4;                                       
            }

            //If not already initialized, initialize particle rest density
            if (m_particleRestDensity == 0.0) {
				float sum = 0.0;
				size_t numFluidCells = 0;

				for (size_t x = 0; x < GRID_WIDTH; ++x) {
                    for (size_t y = 0; y < GRID_HEIGHT; ++y) {
                        const Cell<GRID_CELL_SIZE>& cell = GetCell(x,y );
                        if (cell.GetCellState() == CellState::Fluid)
                        {
                            sum += m_densities[x][y];
						    numFluidCells++;
                        }
                    }
				}

				if (numFluidCells > 0)
					m_particleRestDensity = sum / numFluidCells;
			}
        }

        void EnforceBoundaryConditions() {
            for (size_t x = 0; x < GRID_WIDTH; ++x) {
                for (size_t y = 0; y < GRID_HEIGHT; ++y) {
                    const Cell<GRID_CELL_SIZE>& cell = GetCell(x, y);
                    if (cell.GetCellState() == CellState::Solid){
                         continue;
                    }

                    CellVelocities v = GetCellVelocities(x, y);
                    if (!cell.IsLeftSurfaceCompressible())  v.leftVelocity  = 0.0f;
                    if (!cell.IsRightSurfaceCompressible()) v.rightVelocity = 0.0f;
                    if (!cell.IsUpSurfaceCompressible())    v.upVelocity    = 0.0f;
                    if (!cell.IsDownSurfaceCompressible())  v.downVelocity  = 0.0f;
                }
            }
        }


        CellVelocities GetCellVelocities(int16_t xIndex, int16_t yIndex) {
             if(!ValidCell(xIndex, yIndex)){
                throw std::runtime_error("Attempting to access invalid cell");
            }

            float& leftVelocity = m_horizontalVelocities[(yIndex * (GRID_WIDTH + 1)) + xIndex];
            float& rightVelocity = m_horizontalVelocities[(yIndex * (GRID_WIDTH + 1)) + (xIndex + 1)];
            float& upVelocity = m_verticalVelocities[(xIndex * (GRID_HEIGHT + 1)) + yIndex];
            float& downVelocity = m_verticalVelocities[(xIndex * (GRID_HEIGHT + 1)) + (yIndex + 1)];

            return CellVelocities{leftVelocity, rightVelocity, upVelocity, downVelocity};
        }

        void SimulateParticles(float dt, Vector2 gravity){
            for (auto& particle : m_particles){
                particle.Simulate(dt, gravity,
                    Vector2{static_cast<float>(GRID_CELL_SIZE), static_cast<float>(GRID_CELL_SIZE)},
                    Vector2{static_cast<float>((GRID_WIDTH - 1) * GRID_CELL_SIZE), static_cast<float>((GRID_HEIGHT - 1) * GRID_CELL_SIZE)});

                //Push particles out of obstacles
                const Cell<GRID_CELL_SIZE>& cell = GetCellFromGlobalPos(particle.position);
                if (cell.GetCellState() == CellState::Solid){
                    std::pair<Vector2, std::pair<Direction, Direction>> distanceData = cell.DistancesFromClosestBoundaries(particle.position);

                    auto distance = distanceData.first;
                    auto horizontalDir = distanceData.second.first;
                    auto verticalDir = distanceData.second.second;

                    Direction pushDirection = (distance.x < distance.y) ? horizontalDir : verticalDir;
                    Direction alternatePushDirection = (distance.x >= distance.y) ? horizontalDir : verticalDir;

                    Direction actualPushDirection = pushDirection;
                    if (cell.IsSurfaceCompresssible(pushDirection)){ actualPushDirection = pushDirection;}
                    else if (cell.IsSurfaceCompresssible(alternatePushDirection)){ actualPushDirection = alternatePushDirection;}
                    else {
                        // Emergency fallback: clamp to fluid domain and kill velocity
                        static constexpr float minBound = static_cast<float>(GRID_CELL_SIZE);
                        static constexpr float maxBoundX = static_cast<float>((GRID_WIDTH - 1) * GRID_CELL_SIZE);
                        static constexpr float maxBoundY = static_cast<float>((GRID_HEIGHT - 1) * GRID_CELL_SIZE);
                        particle.position.x = std::clamp(particle.position.x, minBound, maxBoundX);
                        particle.position.y = std::clamp(particle.position.y, minBound, maxBoundY);
                        particle.velocity = {0.0f, 0.0f};
                        continue;
                    }

                    switch (actualPushDirection) {
                        case Direction::LEFT:
                            particle.Move({-distance.x, 0.0f});
                            break;
                        case Direction::RIGHT:
                            particle.Move({distance.x, 0.0f});
                            break;
                        case Direction::UP:
                            particle.Move({0.0f, -distance.y});
                            break;
                        case Direction::DOWN:
                            particle.Move({0.0f, distance.y});
                            break;
                    }                               
                }
            }
        }

        void MakeIncompressible(float dt){
            // Optimal SOR parameter: ω = 2 / (1 + sin(π/N)) where N is
            // the largest grid dimension. Iteration count ~2N for full
            // convergence of the hydrostatic pressure gradient.
            static constexpr size_t N = GRID_WIDTH > GRID_HEIGHT ? GRID_WIDTH : GRID_HEIGHT;
            static constexpr float PI = 3.14159265358979f;
            static constexpr float OVERRELAXATION_CONSTANT = 2.0f / (1.0f + std::sin(PI / N));
            static constexpr size_t DIVERGENCE_ITERATIONS = N * 2;
            m_pressures = {{}};

            for(int iteration = 0; iteration < DIVERGENCE_ITERATIONS; ++iteration){
                //Alternate sweeps
                bool reverse = (iteration % 2) == 1;
                // All the below nonsense is doing is inverting the indexes if we're on a reverse sweep
                for (size_t ii = 0; ii < GRID_WIDTH; ++ii){
                    const size_t x = reverse ? (GRID_WIDTH - 1 - ii) : ii;
                    for(size_t jj = 0; jj < GRID_HEIGHT; jj++){
                        const size_t y = reverse ? (GRID_HEIGHT - 1 - jj) : jj;

                        /* Only fluid cells participate in the pressure solve */
                        const Cell<GRID_CELL_SIZE>& cell = GetCell(x, y);
                        if (cell.GetCellState() != CellState::Fluid){
                            continue;
                        }
                        if (cell.NumCompressibleSurfaces() == 0){
                            continue;
                        }

                        CellVelocities velocities = GetCellVelocities(x, y);

                        const float STIFFNESS = 1.0f;
                        const float divergence = velocities.Divergence() - (STIFFNESS * (m_densities[x][y] - m_particleRestDensity));
                        const float divergenceDelta = (divergence * OVERRELAXATION_CONSTANT)  / cell.NumCompressibleSurfaces();
                        if (cell.IsLeftSurfaceCompressible()) {velocities.leftVelocity += divergenceDelta;}
                        if (cell.IsRightSurfaceCompressible()) {velocities.rightVelocity -= divergenceDelta;}
                        if (cell.IsUpSurfaceCompressible()) {velocities.upVelocity += divergenceDelta;}
                        if (cell.IsDownSurfaceCompressible()) {velocities.downVelocity -= divergenceDelta;}

                        static constexpr float WATER_DENSITY = 1000;
                        m_pressures[x][y] -= divergenceDelta * ((WATER_DENSITY * GRID_CELL_SIZE) / dt);
                    }
                }
            }
        }

            // First is horizontal, second is vertical
        std::pair<StaggeredCornerVelocities, StaggeredCornerVelocities> GridTransferVelocities(Vector2 globalPos ) {
            // For this one, we need to cross grid boundaries, because we are on the staggered grid.
            // Get the cell we're in.
            const Cell<GRID_CELL_SIZE>& cell = GetCellFromGlobalPos(globalPos);

            // Work out which cell we're going to need to be reaching into, it's the closest
            std::pair<Vector2, std::pair<Direction, Direction>> directionData = cell.DistancesFromClosestBoundaries(globalPos);

            // Transfer horizontal velocities
            Direction verticalDir = directionData.second.second;
            // To transfer horizontal velocities, we're reaching into a vertical cell to get the velocity corners
            const Cell<GRID_CELL_SIZE>& verticalCell = cell.GetAdjacentCell(verticalDir);
        
            CellVelocities thisCellVelocities = GetCellVelocities(cell.GetXIndex(), cell.GetYIndex());
            CellVelocities horizontalAdjacentCellVelocities = GetCellVelocities(verticalCell.GetXIndex(), verticalCell.GetYIndex());

            // Deltas must be measured from corner1's position, not from the nearest boundary.
            // corner1 is this cell's left face, at (cellLeft, cellCenterY).
            // deltaX: always from left edge. deltaY: from cell center toward the adjacent cell.
            const float localX = globalPos.x - cell.GetXIndex() * static_cast<float>(GRID_CELL_SIZE);
            const float localY = globalPos.y - cell.GetYIndex() * static_cast<float>(GRID_CELL_SIZE);
            const float deltaXHorizontalVels = localX;
            const float deltaYHorizontalVels = (verticalDir == Direction::DOWN)
                ? (localY - static_cast<float>(GRID_CELL_SIZE) / 2)
                : (static_cast<float>(GRID_CELL_SIZE) / 2 - localY);

            size_t horizontalCorner1FlatIndex = cell.GetYIndex() * (GRID_WIDTH + 1) + cell.GetXIndex();
            size_t horizontalCorner2FlatIndex = cell.GetYIndex() * (GRID_WIDTH + 1) + cell.GetXIndex() + 1;
            size_t horizontalCorner3FlatIndex = verticalCell.GetYIndex() * (GRID_WIDTH + 1) + cell.GetXIndex() + 1;
            size_t horizontalCorner4FlatIndex = verticalCell.GetYIndex() * (GRID_WIDTH + 1) + cell.GetXIndex();

            StaggeredCornerVelocities horizontalCorners = (verticalCell.GetCellState() == CellState::Fluid) ? 
                StaggeredCornerVelocities{
                    &thisCellVelocities.leftVelocity,
                    horizontalCorner1FlatIndex,
                    &thisCellVelocities.rightVelocity,
                    horizontalCorner2FlatIndex,
                    &horizontalAdjacentCellVelocities.rightVelocity,
                    horizontalCorner3FlatIndex,
                    &horizontalAdjacentCellVelocities.leftVelocity,
                    horizontalCorner4FlatIndex,
                    deltaXHorizontalVels,
                    deltaYHorizontalVels} :
                StaggeredCornerVelocities{
                    &thisCellVelocities.leftVelocity,
                    horizontalCorner1FlatIndex,
                    &thisCellVelocities.rightVelocity,
                    horizontalCorner2FlatIndex,
                    nullptr, 
                    std::numeric_limits<size_t>::max(),
                    nullptr, 
                    std::numeric_limits<size_t>::max(),
                    deltaXHorizontalVels, 
                    deltaYHorizontalVels}; 

            // Transfer vertical velocities
            Direction horizontalDir = directionData.second.first;
            // To transfer vertical velocities, we're reaching into a horizontal cell to get the velocity corners
            const Cell<GRID_CELL_SIZE>& horizontalCell = cell.GetAdjacentCell(horizontalDir);
        
            CellVelocities VertcialAdjacentCellVelocities = GetCellVelocities(horizontalCell.GetXIndex(), horizontalCell.GetYIndex());

            // corner1 is this cell's up face, at (cellCenterX, cellTop).
            // corner1→corner2 goes top→bottom (physical Y), so weight formula's dx = localY.
            // corner1→corner4 goes toward adjacent cell (physical X), so weight formula's dy = distance toward adjacent.
            const float deltaXVerticalVels = localY;
            const float deltaYVerticalVels = (horizontalDir == Direction::RIGHT)
                ? (localX - static_cast<float>(GRID_CELL_SIZE) / 2)
                : (static_cast<float>(GRID_CELL_SIZE) / 2 - localX);

            size_t verticalCorner1FlatIndex = cell.GetXIndex() * (GRID_HEIGHT + 1) + cell.GetYIndex();
            size_t verticalCorner2FlatIndex = cell.GetXIndex() * (GRID_HEIGHT + 1) + cell.GetYIndex() + 1;
            size_t verticalCorner3FlatIndex = horizontalCell.GetXIndex() * (GRID_HEIGHT + 1) + cell.GetYIndex() + 1;
            size_t verticalCorner4FlatIndex = horizontalCell.GetXIndex() * (GRID_HEIGHT + 1) + cell.GetYIndex();

            StaggeredCornerVelocities verticalCorners = (horizontalCell.GetCellState() == CellState::Fluid) ? 
                StaggeredCornerVelocities{
                    &thisCellVelocities.upVelocity, 
                    verticalCorner1FlatIndex,
                    &thisCellVelocities.downVelocity, 
                    verticalCorner2FlatIndex,
                    &VertcialAdjacentCellVelocities.downVelocity, 
                    verticalCorner3FlatIndex,
                    &VertcialAdjacentCellVelocities.upVelocity, 
                    verticalCorner4FlatIndex,
                    deltaXVerticalVels, 
                    deltaYVerticalVels} :
                StaggeredCornerVelocities{
                    &thisCellVelocities.upVelocity, 
                    verticalCorner1FlatIndex,
                    &thisCellVelocities.downVelocity, 
                    verticalCorner2FlatIndex,
                    nullptr,
                    std::numeric_limits<size_t>::max(),
                    nullptr, 
                    std::numeric_limits<size_t>::max(),
                    deltaXVerticalVels, 
                    deltaYVerticalVels}; 

            return std::make_pair(horizontalCorners, verticalCorners);
        }

        void VelocityTransferParticlesToGrid(){

            /*
             * For every particle, we work out the weights to the four closest velocity corners on the staggered
             * grid, and we do this for both horizontal and vertical "cells". Then, using the weights, we add the weight*velocities to the
             * staggered grid velocities, whilst also accumulating a counter that just accumulates the weight values for each face. Then
             * right at the end, we perform normalization on each face (ie, each velocity), by just doing a velocity/=accumulatedWeight
             */ 

            /* Zero out current grid velocities, since we're getting them from particles, and weirdnesses with undefined corners can lead to some stray accumulations if you don't do this*/
            m_verticalVelocities = {};
            m_horizontalVelocities = {};

            /* First is the weighted sum (the one that weights velocities), second is just the weight accumulator (just weights, no velocity, used to normalize) */
            std::array<std::pair<float,float>, (GRID_WIDTH + 1) * GRID_HEIGHT> horizontalVelocitiesTransferAccumulations = {};
            std::array<std::pair<float,float>, GRID_WIDTH * (GRID_HEIGHT + 1)> verticalVelocitiesTransferAccumulations = {};
            
            //Accumulate weighted velocities from particles
            for (const auto& particle : m_particles){
                // Get the corners of the staggered "cell" we're transferring velocities into
                auto [horizontalStaggeredCellVelocities, verticalStaggeredCellVelocities] = GridTransferVelocities(particle.position);

                //Compute weights
                const float horizontalWeight1 = (1 - (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (1 - (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float horizontalWeight2 = (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (1 - (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float horizontalWeight3 = (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);
                const float horizontalWeight4 = (1 - (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);

                const float verticalWeight1 = (1 - (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (1 - (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float verticalWeight2 = (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (1 - (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float verticalWeight3 = (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);
                const float verticalWeight4 = (1 - (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);
                
                //Remember that some of these corners may be undefined, if they're adjacent to solid cells or whatnot.
                //It's only corners 3 and 4 that might suffer, since these are the ones that reach into adjacent cells

                //Horizontal velocity accumulation
                horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner1FlatIndex].first += horizontalWeight1 * particle.velocity.x;
                horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner2FlatIndex].first += horizontalWeight2 * particle.velocity.x;
                if (horizontalStaggeredCellVelocities.corner3 != nullptr){
                    horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner3FlatIndex].first += horizontalWeight3 * particle.velocity.x;
                }
                if (horizontalStaggeredCellVelocities.corner4 != nullptr){
                    horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner4FlatIndex].first += horizontalWeight4 * particle.velocity.x;
                }
                
                //Horizontal weight accumulation
                horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner1FlatIndex].second += horizontalWeight1;
                horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner2FlatIndex].second += horizontalWeight2;
                if (horizontalStaggeredCellVelocities.corner3 != nullptr){
                    horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner3FlatIndex].second += horizontalWeight3;
                }
                if (horizontalStaggeredCellVelocities.corner4 != nullptr){
                    horizontalVelocitiesTransferAccumulations[horizontalStaggeredCellVelocities.corner4FlatIndex].second += horizontalWeight4;
                }

                //Vertical velocity accumulation
                verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner1FlatIndex].first += verticalWeight1 * particle.velocity.y;
                verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner2FlatIndex].first += verticalWeight2 * particle.velocity.y;
                if (verticalStaggeredCellVelocities.corner3 != nullptr){
                    verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner3FlatIndex].first += verticalWeight3 * particle.velocity.y;
                }
                if (verticalStaggeredCellVelocities.corner4 != nullptr){
                    verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner4FlatIndex].first += verticalWeight4 * particle.velocity.y;
                }

                //Vertical weight accumulation
                verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner1FlatIndex].second += verticalWeight1;
                verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner2FlatIndex].second += verticalWeight2;
                if (verticalStaggeredCellVelocities.corner3 != nullptr){
                    verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner3FlatIndex].second += verticalWeight3;
                }
                if (verticalStaggeredCellVelocities.corner4 != nullptr){
                    verticalVelocitiesTransferAccumulations[verticalStaggeredCellVelocities.corner4FlatIndex].second += verticalWeight4;
                }
            }

            //Normalize weighted velocities
            for (size_t i = 0; i < m_horizontalVelocities.size(); ++i) {
                if (horizontalVelocitiesTransferAccumulations[i].second > 0.0f) {
                    m_horizontalVelocities[i] = horizontalVelocitiesTransferAccumulations [i].first / horizontalVelocitiesTransferAccumulations[i].second;
                }
            }
            for (size_t i = 0; i < m_verticalVelocities.size(); ++i) {
                if (verticalVelocitiesTransferAccumulations[i].second > 0.0f) {
                    m_verticalVelocities[i] = verticalVelocitiesTransferAccumulations[i].first / verticalVelocitiesTransferAccumulations[i].second;
                }
            }
        }

        float CalcWeightedSumNormalizedParticleVelocityComponent(float corner1, float corner2, float corner3, float corner4,  float weight1, float weight2, float weight3, float weight4, bool corner3Defined, bool corner4Defined){
            const float component1 = weight1 * corner1;
            const float component2 = weight2 * corner2;
            const float component3 = corner3Defined ? weight3 * corner3 : 0.0f;
            const float component4 = corner4Defined ? weight4 * corner4 : 0.0f;

            return (component1 + component2 + component3 + component4) / 
                    (weight1 + weight2 + (corner3Defined ? weight3 : 0.0f) + (corner4Defined ? weight4 : 0.0f));
        }

        void VelocityTransferGridtoParticles(const std::array<float, (GRID_WIDTH + 1) * GRID_HEIGHT>& oldHorizontalVelocities, const std::array<float, GRID_WIDTH * (GRID_HEIGHT + 1)>& oldVerticalVelocities){
            for (auto& particle : m_particles){
                auto [horizontalStaggeredCellVelocities, verticalStaggeredCellVelocities] = GridTransferVelocities(particle.position);

                //Compute weights
                const float horizontalWeight1 = (1 - (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (1 - (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float horizontalWeight2 = (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (1 - (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float horizontalWeight3 = (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);
                const float horizontalWeight4 = (1 - (horizontalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (horizontalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);

                const float verticalWeight1 = (1 - (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (1 - (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float verticalWeight2 = (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (1 - (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE));
                const float verticalWeight3 = (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE) * (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);
                const float verticalWeight4 = (1 - (verticalStaggeredCellVelocities.deltaX / GRID_CELL_SIZE)) * (verticalStaggeredCellVelocities.deltaY  / GRID_CELL_SIZE);

                const bool horizontalCorner3Defined = horizontalStaggeredCellVelocities.corner3 != nullptr;
                const bool horizontalCorner4Defined = horizontalStaggeredCellVelocities.corner4 != nullptr;
                const bool verticalCorner3Defined = verticalStaggeredCellVelocities.corner3 != nullptr;
                const bool verticalCorner4Defined = verticalStaggeredCellVelocities.corner4 != nullptr;

                const float interpolatedNewX = CalcWeightedSumNormalizedParticleVelocityComponent(*horizontalStaggeredCellVelocities.corner1, *horizontalStaggeredCellVelocities.corner2, horizontalCorner3Defined ? *horizontalStaggeredCellVelocities.corner3 : 0.0f, horizontalCorner4Defined ? *horizontalStaggeredCellVelocities.corner4 : 0.0f,
                                                                                        horizontalWeight1, horizontalWeight2, horizontalWeight3, horizontalWeight4,
                                                                                        horizontalCorner3Defined, horizontalCorner4Defined);
                
                const float interpolatedNewY = CalcWeightedSumNormalizedParticleVelocityComponent(*verticalStaggeredCellVelocities.corner1, *verticalStaggeredCellVelocities.corner2, verticalCorner3Defined ? *verticalStaggeredCellVelocities.corner3 : 0.0f, verticalCorner4Defined ? *verticalStaggeredCellVelocities.corner4 : 0.0f,
                                                                                        verticalWeight1, verticalWeight2, verticalWeight3, verticalWeight4,
                                                                                        verticalCorner3Defined, verticalCorner4Defined);

                //Do the same interpolation calculation but with the old, pre incompressible solve velocities, such that we can add the changes to each particle, reducing viscocity (FLIP)
                const float oldHorizontalCorner1 = oldHorizontalVelocities[horizontalStaggeredCellVelocities.corner1FlatIndex];
                const float oldHorizontalCorner2 = oldHorizontalVelocities[horizontalStaggeredCellVelocities.corner2FlatIndex];
                const float oldHorizontalCorner3 = horizontalCorner3Defined ? oldHorizontalVelocities[horizontalStaggeredCellVelocities.corner3FlatIndex] : 0.0f;
                const float oldHorizontalCorner4 = horizontalCorner4Defined ? oldHorizontalVelocities[horizontalStaggeredCellVelocities.corner4FlatIndex] : 0.0f;
                
                const float oldVerticalCorner1 = oldVerticalVelocities[verticalStaggeredCellVelocities.corner1FlatIndex];
                const float oldVerticalCorner2 = oldVerticalVelocities[verticalStaggeredCellVelocities.corner2FlatIndex];
                const float oldVerticalCorner3 = verticalCorner3Defined ? oldVerticalVelocities[verticalStaggeredCellVelocities.corner3FlatIndex] : 0.0f;
                const float oldVerticalCorner4 = verticalCorner4Defined ? oldVerticalVelocities[verticalStaggeredCellVelocities.corner4FlatIndex] : 0.0f;

                const float interpolatedOldX = CalcWeightedSumNormalizedParticleVelocityComponent(oldHorizontalCorner1, oldHorizontalCorner2, oldHorizontalCorner3, oldHorizontalCorner4,
                                                                                                  horizontalWeight1, horizontalWeight2, horizontalWeight3, horizontalWeight4,
                                                                                                  horizontalCorner3Defined, horizontalCorner4Defined);

                const float interpolatedOldY = CalcWeightedSumNormalizedParticleVelocityComponent(oldVerticalCorner1, oldVerticalCorner2, oldVerticalCorner3, oldVerticalCorner4,
                                                                                                  verticalWeight1, verticalWeight2, verticalWeight3, verticalWeight4,
                                                                                                  verticalCorner3Defined, verticalCorner4Defined);

                static constexpr float FLIP_RATIO = 0.75f;
                const float flipVx = particle.velocity.x + (interpolatedNewX - interpolatedOldX);
                const float flipVy = particle.velocity.y + (interpolatedNewY - interpolatedOldY);
                particle.velocity.x = FLIP_RATIO * flipVx + (1.0f - FLIP_RATIO) * interpolatedNewX;
                particle.velocity.y = FLIP_RATIO * flipVy + (1.0f - FLIP_RATIO) * interpolatedNewY;

                // Clamp velocity components pointing into solid walls
                const Cell<GRID_CELL_SIZE>& cell = GetCellFromGlobalPos(particle.position);
                if (!cell.IsLeftSurfaceCompressible() && particle.velocity.x < 0.0f)  particle.velocity.x = 0.0f;
                if (!cell.IsRightSurfaceCompressible() && particle.velocity.x > 0.0f) particle.velocity.x = 0.0f;
                if (!cell.IsUpSurfaceCompressible() && particle.velocity.y < 0.0f)    particle.velocity.y = 0.0f;
                if (!cell.IsDownSurfaceCompressible() && particle.velocity.y > 0.0f)  particle.velocity.y = 0.0f;
            }
        }


        void InitCells(){
            for(size_t x = 0; x < GRID_WIDTH; ++x){
                for(size_t y = 0; y < GRID_HEIGHT; ++y){
                    m_cells[x][y] = std::make_unique<Cell<GRID_CELL_SIZE>>(x, y);
                }
            }
        }

        void HookupCellConnections(){
            for(int16_t x = 0; x < GRID_WIDTH; ++x){
                for(int16_t y = 0; y < GRID_HEIGHT; ++y){

                    Cell<GRID_CELL_SIZE>* leftCell = ValidCell(x - 1, y) ? m_cells[x - 1][y].get() : nullptr;
                    Cell<GRID_CELL_SIZE>* rightCell = ValidCell(x + 1, y) ? m_cells[x + 1][y].get() : nullptr;
                    Cell<GRID_CELL_SIZE>* upCell = ValidCell(x, y - 1) ? m_cells[x][y - 1].get() : nullptr;
                    Cell<GRID_CELL_SIZE>* downCell = ValidCell(x, y + 1) ? m_cells[x][y + 1].get() : nullptr;
                    m_cells[x][y]->HookupAdjacents(leftCell, rightCell, upCell, downCell);
                }
            }
        }

        //The edges of the domain are solid
        void SetSolidBorders(){
            for(int16_t x = 0; x < GRID_WIDTH; ++x){
                m_cells[x][0]->SetCellState(CellState::Solid);
                m_cells[x][GRID_HEIGHT - 1]->SetCellState(CellState::Solid);
            }

            for(int16_t y = 0; y < GRID_HEIGHT; ++y){
                m_cells[0][y]->SetCellState(CellState::Solid);
                m_cells[GRID_WIDTH - 1][y]->SetCellState(CellState::Solid);
            }
        }

        void SeedParticles(){
            for(int16_t x = 0; x < GRID_WIDTH; ++x){
                for(int16_t y = 0; y < GRID_HEIGHT; ++y){
                    if (GetCell(x, y).GetCellState() == CellState::Solid) continue;
                    if (y < static_cast<int16_t>(GRID_HEIGHT / 2.0f)) continue;

                    for(size_t p = 0; p < PARTICLES_PER_CELL; ++p){
                        float jitterX = static_cast<float>(p % 2) * 0.5f + 0.25f;
                        float jitterY = static_cast<float>(p / 2) * 0.5f + 0.25f;
                        if (PARTICLES_PER_CELL == 3) {
                            // Spread 3 particles in a triangle-ish pattern
                            constexpr float offsets[3][2] = {{0.25f, 0.25f}, {0.75f, 0.25f}, {0.5f, 0.75f}};
                            jitterX = offsets[p][0];
                            jitterY = offsets[p][1];
                        }

                        m_particles.emplace_back(
                            Vector2{(x * static_cast<float>(GRID_CELL_SIZE)) + jitterX * GRID_CELL_SIZE,
                                    (y * static_cast<float>(GRID_CELL_SIZE)) + jitterY * GRID_CELL_SIZE},
                            Vector2{0.0f, 0.0f});
                    }
                }
            }
        }
        

        /* Stored in the centers of the vertical cell faces (left/right) */
        std::array<float, (GRID_WIDTH + 1) * GRID_HEIGHT> m_horizontalVelocities = {};
        /* Stored in the centers of the horizontlal cell faces (top/bottom) */
        std::array<float, GRID_WIDTH * (GRID_HEIGHT + 1)> m_verticalVelocities = {};
        /* The +1's are because we're storing on the edge, so we need both the left/right and top/bottom "caps", so to speak */

        std::array<std::array<std::unique_ptr<Cell<GRID_CELL_SIZE>>, GRID_HEIGHT>, GRID_WIDTH> m_cells = {};
        std::array<std::array<float, GRID_HEIGHT>, GRID_WIDTH> m_pressures = {};
        std::array<std::array<float, GRID_HEIGHT>, GRID_WIDTH> m_densities = {};
        float m_particleRestDensity = 0.0f;

        // List of cell-adjacent particles we build each timestep, to do push-apart
        std::array<std::array<std::vector<Particle*>, GRID_HEIGHT>, GRID_WIDTH> m_particlesInCells = {};

        static constexpr size_t PARTICLES_PER_CELL = 4;
        std::vector<Particle> m_particles = {};
        size_t m_debugFrames = 0;
    };



} // namespace sim
