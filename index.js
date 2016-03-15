const sm_a = 6378137;
const sm_b = 6356752.314;

const utm_scale_factor = 0.9996;

const n = (sm_a - sm_b) / (sm_a + sm_b);
const alpha = ((sm_a + sm_b) / 2) * (1 + (n**2 / 4) + (n**4 / 64));
const beta = (-3 * n / 2) + (9 * n**3 / 16) + (-3 * n**5 / 32);
const gamma = (15 * n**2 / 16) + (-15 * n**4 / 32);
const delta = (-35 * n**3 / 48) + (105 * n**5 / 256);
const epsilon = (315 * n**4 / 512);

const alpha_ = ((sm_a + sm_b) / 2) * (1 + (n**2 / 4) + (n**4 / 64));
const beta_ = (3 * n / 2) + (-27 * n**3 / 32) + (269 * n**5 / 512);
const gamma_ = (21 * n**2 / 16) + (-55 * n**4 / 32);
const delta_ = (151 * n**3 / 96) + (-417 * n**5 / 128);
const epsilon_ = (1097 * n**4 / 512);

function ArcLengthOfMeridian(phi) {
    return alpha
        * (phi + (beta * Math.sin(2 * phi))
        + (gamma * Math.sin(4 * phi))
        + (delta * Math.sin(6 * phi))
        + (epsilon * Math.sin(8 * phi)));
}

function FootpointLatitude(y) {
    const y_ = y / alpha_;
    return y_ + (beta_ * Math.sin (2 * y_))
        + (gamma_ * Math.sin (4 * y_))
        + (delta_ * Math.sin (6 * y_))
        + (epsilon_ * Math.sin (8 * y_));
}

function GPS2UTM(phi, lambda, zone) {
    phi = phi / 180 * Math.PI;
    lambda = lambda / 180 * Math.PI;
    const lambda0 = (-183 + (zone * 6)) / 180 * Math.PI;

    const ep2 = (sm_a**2 - sm_b**2) / sm_b**2;
    const nu2 = ep2 * Math.cos(phi)**2;
    const N = sm_a**2 / (sm_b * Math.sqrt(1 + nu2));
    const t = Math.tan (phi);
    const t2 = t * t;
    const l = lambda - lambda0;

    const l3coef = 1 - t2 + nu2;
    const l4coef = 5 - t2 + 9 * nu2 + 4 * (nu2 * nu2);
    const l5coef = 5 - 18 * t2 + (t2 * t2) + 14 * nu2 - 58 * t2 * nu2;
    const l6coef = 61 - 58 * t2 + (t2 * t2) + 270 * nu2 - 330 * t2 * nu2;
    const l7coef = 61 - 479 * t2 + 179 * (t2 * t2) - (t2 * t2 * t2);
    const l8coef = 1385 - 3111 * t2 + 543 * (t2 * t2) - (t2 * t2 * t2);

    let x = N * Math.cos(phi) * l
            + (N / 6 * Math.cos(phi)**3 * l3coef * l**3)
            + (N / 120 * Math.cos(phi)**5 * l5coef * l**5)
            + (N/5040 * Math.cos(phi)**7 * l7coef * l**7);

    let y = ArcLengthOfMeridian(phi)
            + (t / 2 * N * Math.cos (phi)**2 * l**2)
            + (t / 24 * N * Math.cos (phi)**4 * l4coef * l**4)
            + (t / 720 * N * Math.cos (phi)**6 * l6coef * l**6)
            + (t / 40320 * N * Math.cos (phi)**8 * l8coef * l**8);

    x = x * utm_scale_factor + 500000;
    y = y * utm_scale_factor;
    if (y < 0) {
        y = y + 10000000;
    }

    return [x, y, zone];
}

function UTM2GPS(x, y, zone, south) {
    x -= 500000;
    x /= utm_scale_factor;

    if (south) {
        y -= 10000000;
    }
    y /= utm_scale_factor;

    const lambda0 = (-183 + (zone * 6)) / 180 * Math.PI;

    const phif = FootpointLatitude (y);
    const ep2 = (sm_a**2 - sm_b**2) / sm_b**2;
    const cf = Math.cos(phif);
    const nuf2 = ep2 * cf**2;
    const Nf = sm_a**2 / (sm_b * Math.sqrt (1 + nuf2));
    let Nfpow = Nf;

    const tf = Math.tan(phif);
    const tf2 = tf * tf;
    const tf4 = tf2 * tf2;

    const x1frac = 1 / (Nfpow * cf);

    Nfpow *= Nf;
    const x2frac = tf / (2 * Nfpow);

    Nfpow *= Nf;
    const x3frac = 1 / (6 * Nfpow * cf);

    Nfpow *= Nf;
    const x4frac = tf / (24 * Nfpow);

    Nfpow *= Nf;
    const x5frac = 1 / (120 * Nfpow * cf);

    Nfpow *= Nf;
    const x6frac = tf / (720 * Nfpow);

    Nfpow *= Nf;
    const x7frac = 1 / (5040 * Nfpow * cf);

    Nfpow *= Nf;
    const x8frac = tf / (40320 * Nfpow);

    const x2poly = -1 - nuf2;
    const x3poly = -1 - 2 * tf2 - nuf2;
    const x4poly = 5 + 3 * tf2 + 6 * nuf2 - 6 * tf2 * nuf2 - 3 * (nuf2 *nuf2)
        - 9 * tf2 * (nuf2 * nuf2);
    const x5poly = 5 + 28 * tf2 + 24 * tf4 + 6 * nuf2 + 8 * tf2 * nuf2;
    const x6poly = -61 - 90 * tf2 - 45 * tf4 - 107 * nuf2 + 162 * tf2 * nuf2;
    const x7poly = -61 - 662 * tf2 - 1320 * tf4 - 720 * (tf4 * tf2);
    const x8poly = 1385 + 3633 * tf2 + 4095 * tf4 + 1575 * (tf4 * tf2);

    const phi = (phif + x2frac * x2poly * (x * x)
        + x4frac * x4poly * x**4
        + x6frac * x6poly * x**6
        + x8frac * x8poly * x**8) * 180 / Math.PI;

    const lambda = (lambda0 + x1frac * x
        + x3frac * x3poly * x**3
        + x5frac * x5poly * x**5
        + x7frac * x7poly * x**7) * 180 / Math.PI;
    return [phi, lambda];
}

export {GPS2UTM, UTM2GPS};
