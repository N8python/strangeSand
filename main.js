import Stats from "./stats.js";
const canvas = document.getElementById('canvas');
const ctx = canvas.getContext('2d');

canvas.width = window.innerWidth;
canvas.height = window.innerHeight;
const canvasWidth = canvas.width;
const canvasHeight = canvas.height;

const numSpheres = 12000;
const gravity = 2000; // Gravitational acceleration
const restitution = 0.75; // Coefficient of restitution
const dragCoefficient = 0.1;
const coefficientOfKinematicFriction = 0.01;

const sphereFields = 4; // x, y, vx, vy, radius, mass, oldX, oldY
const X = 0;
const Y = 1;
const OLD_X = 2;
const OLD_Y = 3;
const spheres = new Float32Array(numSpheres * sphereFields);
const sphereRadius = 2.5;
const cellSize = 2 * sphereRadius;
const sumOfRadiiSquared = cellSize * cellSize;
const invCellSize = 1 / cellSize;
const maxParticlesPerCell = 4;
const gridWidth = Math.ceil(canvasWidth / cellSize);
const gridHeight = Math.ceil(canvasHeight / cellSize);
const grid = new Uint32Array(gridWidth * gridHeight * maxParticlesPerCell);
const nextInCell = new Uint8Array(gridWidth * gridHeight);

function resolveCollision(sphere1, sphere2) {
    const sphere1Accessor = sphere1 * sphereFields;
    const sphere2Accessor = sphere2 * sphereFields;
    const sphere1X = spheres[sphere1Accessor + X];
    const sphere1Y = spheres[sphere1Accessor + Y];
    const sphere2X = spheres[sphere2Accessor + X];
    const sphere2Y = spheres[sphere2Accessor + Y];
    const xDist = sphere2X - sphere1X;
    const yDist = sphere2Y - sphere1Y;
    let distance = (xDist * xDist + yDist * yDist);
    if (distance < sumOfRadiiSquared) {
        distance = Math.sqrt(distance);
        const overlap = 0.5 * (cellSize - distance);
        const normalX = overlap * xDist / distance;
        const normalY = overlap * yDist / distance;
        spheres[sphere1Accessor + X] = sphere1X - normalX;
        spheres[sphere1Accessor + Y] = sphere1Y - normalY;
        spheres[sphere2Accessor + X] = sphere2X + normalX;
        spheres[sphere2Accessor + Y] = sphere2Y + normalY;
    } else {
        // Move the spheres towards each other
        /* distance = Math.sqrt(distance);
         const overlap = 0.5 * (cellSize - distance);
         const normalX = overlap * xDist / distance;
         const normalY = overlap * yDist / distance;
         distance += 0.1;
         spheres[sphere1Accessor + X] = sphere1X - (0.1 / distance) * normalX;
         spheres[sphere1Accessor + Y] = sphere1Y - (0.1 / distance) * normalY;
         spheres[sphere2Accessor + X] = sphere2X + (0.1 / distance) * normalX;
         spheres[sphere2Accessor + Y] = sphere2Y + (0.1 / distance) * normalY;*/


    }
}
class HilbertCurve {
    constructor(n) {
        this.n = n;
        this.size = 1 << n; // 2^n, the size of each dimension
        this.totalElements = this.size * this.size; // Total elements in the curve
        this.cache = new Int32Array(this.totalElements);
        this.precompute();
    }

    precompute() {
        for (let y = 0; y < this.size; y++) {
            for (let x = 0; x < this.size; x++) {
                const index = x * this.size + y;
                this.cache[index] = this.computeHilbertIndex(x, y);
            }
        }
    }

    computeHilbertIndex(x, y) {
        let result = 0;
        let s = this.size >> 1;
        while (s > 0) {
            const rx = (x & s) > 0;
            const ry = (y & s) > 0;
            result += s * s * ((3 * rx) ^ ry);
            [x, y] = this.rotateAndFlip(s, x, y, rx, ry);
            s >>= 1;
        }
        return result;
    }

    rotateAndFlip(s, x, y, rx, ry) {
        if (ry === 0) {
            if (rx === 1) {
                x = s - 1 - x;
                y = s - 1 - y;
            }
            [x, y] = [y, x];
        }
        return [x, y];
    }

    getIndex(x, y) {
        // Convert x, y to integers in the range [0, 2^n - 1] and wrap around using modulus
        let xi = Math.floor(x) % this.size;
        let yi = Math.floor(y) % this.size;

        // Calculate cache index
        const cacheIndex = xi * this.size + yi;
        return this.cache[cacheIndex];
    }
}

const hilbert = new HilbertCurve(10);
let tempSphere = new Float32Array(sphereFields);

function swapSpheres(sphere1, sphere2) {
    const sphere1Accessor = sphere1 * sphereFields;
    const sphere2Accessor = sphere2 * sphereFields;
    for (let i = 0; i < sphereFields; i++) {
        tempSphere[i] = spheres[sphere1Accessor + i];
        spheres[sphere1Accessor + i] = spheres[sphere2Accessor + i];
        spheres[sphere2Accessor + i] = tempSphere[i];
    }
}

function inPlaceMerge(spheres, sphereFields, low, mid, high, hilbert, X, Y) {
    let start = low;
    let start2 = mid;
    while (start < start2 && start2 < high) {
        let index = start * sphereFields;
        let nextIndex = start2 * sphereFields;
        if (hilbert.getIndex(spheres[index + X], spheres[index + Y]) >
            hilbert.getIndex(spheres[nextIndex + X], spheres[nextIndex + Y])) {
            // Rotate the elements between start and start2 right by one position
            let tempX = spheres[nextIndex + X];
            let tempY = spheres[nextIndex + Y];
            for (let i = start2; i > start; i--) {
                spheres[i * sphereFields + X] = spheres[(i - 1) * sphereFields + X];
                spheres[i * sphereFields + Y] = spheres[(i - 1) * sphereFields + Y];
            }
            spheres[start * sphereFields + X] = tempX;
            spheres[start * sphereFields + Y] = tempY;
            start2++;
        }
        start++;
    }
}

function inPlaceMergeSort(spheres, sphereFields, low, high, hilbert, X, Y) {
    if (high - low < 2) return; // No need to sort a single element

    let mid = Math.floor((low + high) / 2);
    inPlaceMergeSort(spheres, sphereFields, low, mid, hilbert, X, Y); // Sort the first half
    inPlaceMergeSort(spheres, sphereFields, mid, high, hilbert, X, Y); // Sort the second half
    inPlaceMerge(spheres, sphereFields, low, mid, high, hilbert, X, Y); // Merge both halves
}

function init() {
    let sphereList = [];
    for (let i = 0; i < numSpheres; i++) {
        const radius = sphereRadius;
        sphereList.push({
            x: radius + Math.random() * (canvasWidth - 2 * radius),
            y: radius + Math.random() * (canvasHeight - 2 * radius),
        });
    }
    //sphereList.sort((a, b) => hilbert.getIndex(a.x, a.y) - hilbert.getIndex(b.x, b.y));
    for (let i = 0; i < numSpheres; i++) {
        const index = i * sphereFields;
        spheres[index + X] = sphereList[i].x;
        spheres[index + Y] = sphereList[i].y;
        spheres[index + OLD_X] = sphereList[i].x;
        spheres[index + OLD_Y] = sphereList[i].y;
    }
    requestAnimationFrame(animate);
}
let time = performance.now();
const STEPS = 4;
const stats = new Stats();
document.body.appendChild(stats.dom);
let mouseX = 0,
    mouseY = 0,
    mouseDown = false;
document.addEventListener('mousemove', (e) => {
    if (mouseDown) {
        mouseX = e.clientX;
        mouseY = e.clientY;
    }
});
document.addEventListener('mousedown', (e) => {
    mouseDown = true;
    mouseX = e.clientX;
    mouseY = e.clientY;
});
document.addEventListener('mouseup', () => {
    mouseDown = false;
});

function animate() {
    stats.update();
    ctx.clearRect(0, 0, canvasWidth, canvasHeight);
    let deltaTime = performance.now() - time;
    time = performance.now();
    deltaTime /= 1000; // Convert to seconds
    deltaTime = Math.min(deltaTime, 0.016); // Clamp to 5ms
    let timeSquared = deltaTime * deltaTime;
    let stepTime = deltaTime / STEPS;

    // Clear grid
    for (let i = 0; i < STEPS; i++) {
        // Resolve collisions using the grid and adjacent cells
        let reverse = i % 2 == 0;
        /* for (let gy = 0; gy < gridHeight; gy++) {
             for (let gx = 0; gx < gridWidth; gx++) {*/
        for (let gy = reverse ? gridHeight - 1 : 0; reverse ? gy >= 0 : gy < gridHeight; reverse ? gy-- : gy++) {
            for (let gx = reverse ? gridWidth - 1 : 0; reverse ? gx >= 0 : gx < gridWidth; reverse ? gx-- : gx++) {

                const baseGridIndex = (gy * gridWidth + gx) * maxParticlesPerCell;
                const particlesInBaseCell = nextInCell[gy * gridWidth + gx];

                // Check collisions within the cell
                for (let idx1 = 0; idx1 < particlesInBaseCell; idx1++) {
                    const sphere1 = grid[baseGridIndex + idx1];

                    // Check collisions with other spheres in the same cell
                    for (let idx2 = idx1 + 1; idx2 < particlesInBaseCell; idx2++) {
                        const sphere2 = grid[baseGridIndex + idx2];
                        resolveCollision(sphere1, sphere2);
                    }

                    // Check collisions with spheres in adjacent cells
                    const startX = Math.max(gx - 1, 0);
                    const endX = Math.min(gx + 1, gridWidth - 1);
                    const startY = Math.max(gy - 1, 0);
                    const endY = Math.min(gy + 1, gridHeight - 1);

                    for (let adjY = startY; adjY <= endY; adjY++) {
                        for (let adjX = startX; adjX <= endX; adjX++) {
                            if (adjX == gx && adjY == gy) {
                                continue; // Skip the base cell as it's already handled
                            }

                            const adjGridIndex = (adjY * gridWidth + adjX) * maxParticlesPerCell;
                            const particlesInAdjCell = nextInCell[adjY * gridWidth + adjX];

                            for (let idxAdj = 0; idxAdj < particlesInAdjCell; idxAdj++) {
                                const sphere2 = grid[adjGridIndex + idxAdj];
                                resolveCollision(sphere1, sphere2);
                            }
                        }
                    }
                }
            }
        }
        grid.fill(0);
        nextInCell.fill(0);
        let timeSquared = stepTime * stepTime;
        /*let swaps = 0;
        for (let i = 1; i < numSpheres; i++) {
            const index = i * sphereFields;
            const currentX = spheres[index + X];
            const currentY = spheres[index + Y];
            const currentIndexValue = hilbert.getIndex(currentX, currentY);

            let j = i - 1;
            while (j >= 0) {
                const prevIndex = j * sphereFields;
                const prevIndexValue = hilbert.getIndex(spheres[prevIndex + X], spheres[prevIndex + Y]);

                if (prevIndexValue > currentIndexValue) {
                    swapSpheres(j + 1, j);
                    swaps++;
                } else {
                    break;
                }
                j--;
            }
        }*/
        // Update grid with new positions
        for (let i = 0; i < numSpheres; i++) {
            const index = i * sphereFields;
            const xDiff = spheres[index + X] - spheres[index + OLD_X];
            const yDiff = spheres[index + Y] - spheres[index + OLD_Y];
            let newX = 2 * spheres[index + X] - spheres[index + OLD_X];
            let newY = 2 * spheres[index + Y] - spheres[index + OLD_Y] + 0.5 * gravity * timeSquared;


            if (newX - sphereRadius <= 0) {
                newX = sphereRadius - xDiff * restitution;
            }
            if (newX + sphereRadius >= canvasWidth) {
                newX = canvasWidth - sphereRadius + xDiff * restitution;
            }
            if (newY - sphereRadius <= 0) {
                newY = sphereRadius - yDiff * restitution;
            }
            if (newY + sphereRadius >= canvasHeight) {
                newY = canvasHeight - sphereRadius + yDiff * restitution;
            }
            // Attract to mouse
            if (mouseDown) {
                let dx = mouseX - spheres[index + X];
                let dy = mouseY - spheres[index + Y];
                let distance = (dx * dx + dy * dy);
                if (distance < 40000) {
                    distance = Math.sqrt(distance);
                    const force = 0.01 * (200 - distance);
                    /*   const angle = Math.atan2(dy, dx);
                       newX += Math.cos(angle) * force * deltaTime;
                       newY += Math.sin(angle) * force * deltaTime;*/
                    dx /= distance;
                    dy /= distance;
                    newX += dx * force * deltaTime;
                    newY += dy * force * deltaTime;

                }
            }

            spheres[index + OLD_X] = spheres[index + X];
            spheres[index + OLD_Y] = spheres[index + Y];
            const lerpFactor = 0.999;
            spheres[index + X] = newX * lerpFactor + spheres[index + X] * (1 - lerpFactor);
            spheres[index + Y] = newY * lerpFactor + spheres[index + Y] * (1 - lerpFactor);

            // Update grid
            const gridX = Math.floor(newX / cellSize);
            const gridY = Math.floor(newY / cellSize);
            const gridIndex = (gridY * gridWidth + gridX) * maxParticlesPerCell;
            const nextSlot = nextInCell[gridY * gridWidth + gridX]++;
            if (nextSlot < maxParticlesPerCell) {
                grid[gridIndex + nextSlot] = i;
            }
        }

    }

    // Draw spheres
    ctx.fillStyle = '#C2B280';
    ctx.beginPath();
    for (let i = 0; i < numSpheres; i++) {
        const index = i * sphereFields;
        ctx.moveTo(spheres[index + X] + sphereRadius, spheres[index + Y]);
        ctx.arc(spheres[index + X], spheres[index + Y], sphereRadius, 0, 2 * Math.PI);
    }
    ctx.fill();
    requestAnimationFrame(animate);
}

init();