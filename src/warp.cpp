/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/warp.h>
#include <nori/vector.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

Vector3f Warp::sampleUniformHemisphere(Sampler *sampler, const Normal3f &pole) {
    // Naive implementation using rejection sampling
    Vector3f v;
    do {
        v.x() = 1.f - 2.f * sampler->next1D();
        v.y() = 1.f - 2.f * sampler->next1D();
        v.z() = 1.f - 2.f * sampler->next1D();
    } while (v.squaredNorm() > 1.f);

    if (v.dot(pole) < 0.f)
        v = -v;
    v /= v.norm();

    return v;
}

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    // by inverse sampling, see report
    Point2f warped;
    float r = sqrt(sample.x());
    float theta = 2 * M_PI * sample.y();
    warped.x() = r * cos(theta);
    warped.y() = r * sin(theta);

    return warped;
}

float Warp::squareToUniformDiskPdf(const Point2f &p) {
    if (p.norm() <= 1) {
        return 1.f / M_PI;
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToUniformSphereCap(const Point2f &sample, float cosThetaMax) {
    // by inverse sampling, see report
    float cosTheta = (1 - sample[0]) + sample[0] * cosThetaMax;
    float sinTheta = sqrt(1.f - cosTheta * cosTheta);
    float phi = sample[1] * 2 * M_PI;

    return Vector3f(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,
                    cosTheta);
}

float Warp::squareToUniformSphereCapPdf(const Vector3f &v, float cosThetaMax) {
    if (abs(v.norm() - 1.f) < Epsilon && cosThetaMax <= v.z() && v.z() <= 1.f) {
        return 1.f / (2.f * M_PI * (1 - cosThetaMax));
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    return squareToUniformSphereCap(sample, -1.f);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    if (abs(v.norm() - 1.f) < Epsilon && -1.f <= v.z() && v.z() <= 1.f) {
        return 1.f / (4.f * M_PI);
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    return squareToUniformSphereCap(sample, 0.f);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (abs(v.norm() - 1.f) < Epsilon && 0.f <= v.z() && v.z() <= 1.f) {
        return 1.f / (2.f * M_PI);
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    // by Malley's method
    float r = sqrt(sample.x());
    float phi = 2 * M_PI * sample.y();
    float x = r * cos(phi);
    float y = r * sin(phi);

    return (Vector3f(x, y, sqrt(fmax(0.f, 1.f - pow(x, 2) - pow(y, 2)))));
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    if (abs(v.norm() - 1.f) < Epsilon && 0.f <= v.z() && v.z() <= 1.f) {
        return v.z() / M_PI;
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    // inverse sampling, see report
    float theta = atan(sqrt(-alpha * alpha * log(1-sample[0])));
    float phi = sample[1] * 2.f * M_PI;

    return Vector3f(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (abs(m.norm() - 1.f) < Epsilon && 0 < m.z() && m.z() <= 1) { // m.z() in denominator can't be zero
        float cosTheta = m(2);
        float tanTheta_2 = (1.f - pow(cosTheta, 2)) / pow(cosTheta, 2);
        float numerator = exp(-tanTheta_2 / (alpha * alpha));
        float denominator = M_PI * alpha * alpha * pow(cosTheta, 3);

        return numerator / denominator;
    } else {
        return 0.f;
    }
}

Vector3f Warp::squareToUniformTriangle(const Point2f &sample) {
    float su1 = sqrtf(sample.x());
    float u = 1.f - su1, v = sample.y() * su1;
    return Vector3f(u,v,1.f-u-v);
}

NORI_NAMESPACE_END
