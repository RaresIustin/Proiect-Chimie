const N = 50;
const iter = 4;
const SCALE = 10;
let t = 0;
let densitatea = 1000;
let fluid;
let offsetX, offsetY;

function setup() {
    createCanvas(N * SCALE, N * SCALE);
    console.log("densitate:" + densitatea);
    fluid = new Fluid(2, 0.00000000001, 0.000001);
    offsetX = (width - N * SCALE) / 2;
    offsetY = (height - N * SCALE) / 2;
}

function draw() {
    if (keyIsPressed) {
        if (key >= '0' && key <= '9') {
            densitatea = densitatea + parseInt(key) - '0';
            console.log("densitate noua " + densitatea);
        }
    }

    background(0); 
    let cx = int(0.5 * width / SCALE);
    let cy = int(0.5 * height / SCALE);
    for (let i = -1; i <= 1; i++) {
        for (let j = -1; j <= 1; j++) {
            fluid.addDensity(cx + i, cy + j, densitatea);
        }
    }
    for (let i = 0; i < 2; i++) {
        let angle = noise(t) * TWO_PI * 2;
        let v = p5.Vector.fromAngle(angle);
        v.mult(0.2);
        t += 0.01;
        fluid.addVelocity(cx, cy, v.x, v.y);
    }

    fluid.step();
    fluid.renderD();
}

function mouseDragged() {
    fluid.addDensity((mouseX-offsetX) / SCALE, (mouseY-offsetY) / SCALE, densitatea);
    let amtX = mouseX - pmouseX;
    let amtY = mouseY - pmouseY;
    fluid.addVelocity((mouseX-offsetX) / SCALE, (mouseY-offsetY) / SCALE, amtX, amtY);
}

class Fluid {
    constructor(dt, diffusion, viscosity) {
        this.size = N;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;

        this.s = new Array(N * N).fill(0);
        this.density = new Array(N * N).fill(0);

        this.Vx = new Array(N * N).fill(0);
        this.Vy = new Array(N * N).fill(0);

        this.Vx0 = new Array(N * N).fill(0);
        this.Vy0 = new Array(N * N).fill(0);
    }

    IX(x, y) {
        x = constrain(x, 0, N - 1);
        y = constrain(y, 0, N - 1);
        return x + (y * N);
    }

    addDensity(x, y, amount) {
        let index = this.IX(x, y);
        this.density[index] += amount;
    }

    addVelocity(x, y, amountX, amountY) {
        let index = this.IX(x, y);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
    }

    step() {
        let N = this.size;
        let visc = this.visc;
        let diff = this.diff;
        let dt = this.dt;
        let Vx = this.Vx;
        let Vy = this.Vy;
        let Vx0 = this.Vx0;
        let Vy0 = this.Vy0;
        let s = this.s;
        let density = this.density;

        this.diffuse(1, Vx0, Vx, visc, dt);
        this.diffuse(2, Vy0, Vy, visc, dt);
        this.project(Vx0, Vy0, Vx, Vy);
        this.advect(1, Vx, Vx0, Vx0, Vy0, dt);
        this.advect(2, Vy, Vy0, Vx0, Vy0, dt);
        this.project(Vx, Vy, Vx0, Vy0);
        this.diffuse(0, s, density, diff, dt);
        this.advect(0, density, s, Vx, Vy, dt);
    }

    renderD() {
        colorMode(HSB, 250);
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                let x = i * SCALE + offsetX;
                let y = j * SCALE + offsetY;
                let d = this.density[this.IX(i, j)];
                fill((d + 50) % 25, 60, d);
                noStroke();
                square(x, y, SCALE);
            }
        }
    }

    renderV() {
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                let x = i * SCALE + offsetX;
                let y = j * SCALE + offsetY;
                let vx = this.Vx[this.IX(i, j)];
                let vy = this.Vy[this.IX(i, j)];
                stroke(255);

                if (!(abs(vx) < 0.1 && abs(vy) <= 0.1)) {
                    line(x, y, x + vx * SCALE, y + vy * SCALE);
                }
            }
        }
    }

    fadeD() {
        for (let i = 0; i < this.density.length; i++) {
            let d = this.density[i];
            this.density[i] = constrain(d - 0.02, 0, 255);
        }
    }

    diffuse(b, x, x0, diff, dt) {
        let a = dt * diff * (N - 2) * (N - 2);
        this.lin_solve(b, x, x0, a, 1 + 4 * a);
    }

    lin_solve(b, x, x0, a, c) {
        let cRecip = 1.0 / c;
        for (let k = 0; k < iter; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    x[this.IX(i, j)] =
                        (x0[this.IX(i, j)] +
                        a * (x[this.IX(i + 1, j)] +
                            x[this.IX(i - 1, j)] +
                            x[this.IX(i, j + 1)] +
                            x[this.IX(i, j - 1)])) * cRecip;
                }
            }
            this.set_bnd(b, x);
        }
    }

    project(velocX, velocY, p, div) {
        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                div[this.IX(i, j)] = -0.5 * (
                    velocX[this.IX(i + 1, j)] -
                    velocX[this.IX(i - 1, j)] +
                    velocY[this.IX(i, j + 1)] -
                    velocY[this.IX(i, j - 1)]
                ) / N;
                p[this.IX(i, j)] = 0;
            }
        }

        this.set_bnd(0, div);
        this.set_bnd(0, p);
        this.lin_solve(0, p, div, 1, 4);

        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                velocX[this.IX(i, j)] -= 0.5 * (p[this.IX(i + 1, j)] - p[this.IX(i - 1, j)]) * N;
                velocY[this.IX(i, j)] -= 0.5 * (p[this.IX(i, j + 1)] - p[this.IX(i, j - 1)]) * N;
            }
        }
        this.set_bnd(1, velocX);
        this.set_bnd(2, velocY);
    }

    advect(b, d, d0, velocX, velocY, dt) {
        let i0, i1, j0, j1;

        let dtx = dt * (N - 2);
        let dty = dt * (N - 2);

        let s0, s1, t0, t1;
        let tmp1, tmp2, x, y;

        let Nfloat = N;
        let ifloat, jfloat;
        let i, j;

        for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
            for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.IX(i, j)];
                tmp2 = dty * velocY[this.IX(i, j)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;
                if (x < 0.5) x = 0.5;
                if (x > Nfloat + 0.5) x = Nfloat + 0.5;
                i0 = floor(x);
                i1 = i0 + 1.0;
                if (y < 0.5) y = 0.5;
                if (y > Nfloat + 0.5) y = Nfloat + 0.5;
                j0 = floor(y);
                j1 = j0 + 1.0;

                s1 = x - i0;
                s0 = 1.0 - s1;
                t1 = y - j0;
                t0 = 1.0 - t1;

                let i0i = int(i0);
                let i1i = int(i1);
                let j0i = int(j0);
                let j1i = int(j1);

                d[this.IX(i, j)] =
                    s0 * (t0 * d0[this.IX(i0i, j0i)] + t1 * d0[this.IX(i0i, j1i)]) +
                    s1 * (t0 * d0[this.IX(i1i, j0i)] + t1 * d0[this.IX(i1i, j1i)]);
            }
        }
        this.set_bnd(b, d);
    }

    set_bnd(b, x) {
        for (let i = 1; i < N - 1; i++) {
            x[this.IX(i, 0)] = b == 2 ? -x[this.IX(i, 1)] : x[this.IX(i, 1)];
            x[this.IX(i, N - 1)] = b == 2 ? -x[this.IX(i, N - 2)] : x[this.IX(i, N - 2)];
        }
        for (let j = 1; j < N - 1; j++) {
            x[this.IX(0, j)] = b == 1 ? -x[this.IX(1, j)] : x[this.IX(1, j)];
            x[this.IX(N - 1, j)] = b == 1 ? -x[this.IX(N - 2, j)] : x[this.IX(N - 2, j)];
        }

        x[this.IX(0, 0)] = 0.5 * (x[this.IX(1, 0)] + x[this.IX(0, 1)]);
        x[this.IX(0, N - 1)] = 0.5 * (x[this.IX(1, N - 1)] + x[this.IX(0, N - 2)]);
        x[this.IX(N - 1, 0)] = 0.5 * (x[this.IX(N - 2, 0)] + x[this.IX(N - 1, 1)]);
        x[this.IX(N - 1, N - 1)] = 0.5 * (x[this.IX(N - 2, N - 1)] + x[this.IX(N - 1, N - 2)]);
    }
}