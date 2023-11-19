/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Romain Prévost

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

#include <nori/shape.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/warp.h>
#include <nori/common.h>

NORI_NAMESPACE_BEGIN

class Sphere : public Shape {
public:
    Sphere(const PropertyList & propList) {
        m_position = propList.getPoint3("center", Point3f());
        m_radius = propList.getFloat("radius", 1.f);

        m_bbox.expandBy(m_position - Vector3f(m_radius));
        m_bbox.expandBy(m_position + Vector3f(m_radius));
    }

    virtual BoundingBox3f getBoundingBox(uint32_t index) const override { return m_bbox; }

    virtual Point3f getCentroid(uint32_t index) const override { return m_position; }

    virtual bool rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const override {
        /* compute a, b and c */
        Vector3f oc = ray.o - m_position;
        float a = ray.d.squaredNorm();
        float b = 2.0 * oc.dot(ray.d);
        float c = oc.squaredNorm() - m_radius * m_radius;

        /* compute and check discriminant */
        float discr = b * b - 4.0 * a * c;
        if (discr < 0) {
            return false;
        }

        /* compute smallest real positive root */
        float sqrtd = sqrt(discr);
        float root1 = (-b + sqrtd) / (2.0 * a);
        float root2 = (-b - sqrtd) / (2.0 * a);
        if (root2 < ray.mint || root2 > ray.maxt) {
            if (root1 < ray.mint || root1 > ray.maxt) {
                return false;
            }
            t = root1;
            return true;
        } else {
            t = root2;
            return true;
        }
    }

    virtual void setHitInformation(uint32_t index, const Ray3f &ray, Intersection & its) const override {
        /* fill its (uv coordinates are computed based on lecture 2) */
        its.p = ray(its.t);
        Point3f pc = (its.p - m_position).normalized();
        its.uv[0] = std::atan2(pc.y(), pc.x());
        its.uv[1] = std::asin(pc.z());
        its.uv[0] = 0.5 + its.uv[0] / (2.f * M_PI);
        its.uv[1] = 0.5 + its.uv[1] / M_PI;
        its.shFrame = Frame((its.p - m_position).normalized());
        its.geoFrame = Frame((its.p - m_position).normalized());
        its.mesh = this;
    }

    virtual void sampleSurface(ShapeQueryRecord & sRec, const Point2f & sample) const override {
        Vector3f q = Warp::squareToUniformSphere(sample);
        sRec.p = m_position + m_radius * q;
        sRec.n = q;
        sRec.pdf = std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }
    virtual float pdfSurface(const ShapeQueryRecord & sRec) const override {
        return std::pow(1.f/m_radius,2) * Warp::squareToUniformSpherePdf(Vector3f(0.0f,0.0f,1.0f));
    }


    virtual std::string toString() const override {
        return tfm::format(
                "Sphere[\n"
                "  center = %s,\n"
                "  radius = %f,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "]",
                m_position.toString(),
                m_radius,
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"));
    }

protected:
    Point3f m_position;
    float m_radius;
};

NORI_REGISTER_CLASS(Sphere, "sphere");
NORI_NAMESPACE_END
