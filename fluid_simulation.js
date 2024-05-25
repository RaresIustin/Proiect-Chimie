import { Fluid } from './Fluid.js';

const N = 200;
const iter = 16;
const SCALE = 4;
let t = 0;
let densitatea = 200;
let fluid;

function setup() {
    createCanvas(N * SCALE, N * SCALE).id('fluid-canvas');
    console.log("densitate:");
    fluid = new Fluid(0.2, 0, 0.0000001);
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
    fluid.addDensity(mouseX / SCALE, mouseY / SCALE, densitatea);
    let amtX = mouseX - pmouseX;
    let amtY = mouseY - pmouseY;
    fluid.addVelocity(mouseX / SCALE, mouseY / SCALE, amtX, amtY);
}
