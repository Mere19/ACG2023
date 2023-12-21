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

#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/scene.h>
#include <nori/camera.h>
#include <nori/photon.h>

NORI_NAMESPACE_BEGIN

/// progressive photon mapper
/// ref: http://graphics.ucsd.edu/~henrik/papers/progressive_photon_mapping/progressive_photon_mapping.pdf
struct HitPoint {
public:
    Point3f x;      // hit location
    Vector3f n;     // normal at x
    Vector3f w;     // ray direction
    const BSDF *bsdf;     // bsdf at x
    Point2i xy;     // pixel location
    Point2f uv;     // uv location
    Color3f t;    // weight
    float R;        // current photon radius
    int N;          // accumulated photon count
    Color3f tao;    // accumulated reflected flux
};

class ProgressivePhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;

    ProgressivePhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.f);
        m_iterCount = props.getInteger("iterCount", 16);
        m_alpha = props.getFloat("alpha", 0.9);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount * m_iterCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon maps */
        for (int i = 0; i < m_iterCount; i ++) {
            m_photonMaps.push_back(std::unique_ptr<PhotonMap>(new PhotonMap()));
            m_photonMaps[i]->reserve(m_photonCount);
        }

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;

        /* multiple photon passes */
        m_emittedPhotonCount = 0;
        for (int i = 0; i < m_iterCount; i ++) {
            /* trace photons */
            int currPhotonCount = 0;
            while (currPhotonCount < m_photonCount) {
                /* sample photon */
                const Emitter* light = scene->getRandomEmitter(sampler->next1D());
                const int n_lights = scene->getLights().size();
                Ray3f currRay;
                Color3f t = 1.f;
                Color3f w = light->samplePhoton(currRay, sampler->next2D(), sampler->next2D()) * n_lights;
                Intersection its;
                m_emittedPhotonCount ++;
                
                /* trace photon */
                while (true) {
                    /* return black when no intersection */
                    if (!scene->rayIntersect(currRay, its))
                        break;

                    /* if diffuse surface, add photon record */
                    if (its.mesh->getBSDF()->isDiffuse()) {
                        m_photonMaps[i]->push_back(Photon(
                            its.p,
                            -currRay.d,
                            t * w
                        ));
                        currPhotonCount ++;
                        if (currPhotonCount == m_photonCount) {
                            break;
                        }
                    }

                    /* Russian roulette */
                    float p = std::min(t.maxCoeff(), 0.99f);
                    if (sampler->next1D() > p || p <= 0.f) {
                        break;
                    }
                    t /= p;

                    /* sample brdf */
                    BSDFQueryRecord bRec(its.toLocal(-currRay.d));
                    bRec.p = its.p;
                    bRec.uv = its.uv;
                    Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                    /* recursion */
                    t *= brdf;
                    currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
                }
            }
            m_emittedPhotonCounts.push_back(m_emittedPhotonCount);

            /* Build the photon map */
            m_photonMaps[i]->build();
        }
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
        Intersection its;
        Ray3f currRay = _ray;
        Color3f t = 1.f;
        Color3f color = 0;
        HitPoint hp;

        /* ray trace */ 
        while (true) {
            if (!scene->rayIntersect(currRay, its)) {
                break;
            }

            /* compute Le when intersection is an emitter */
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord lRec(currRay.o, its.p, its.shFrame.n);
                color += t * its.mesh->getEmitter()->eval(lRec);
            }

            /* store hit point for diffuse surfaces */
            if (its.mesh->getBSDF()->isDiffuse()) {
                hp.x = its.p;
                hp.n = its.shFrame.n;
                hp.w = its.toLocal(-currRay.d);
                hp.bsdf = its.mesh->getBSDF();
                hp.uv = its.uv;
                hp.t = t;
                hp.R = m_photonRadius;
                hp.N = 0;
                hp.tao = 0;

                break;
            }

            /* Russian roulette */
            float p = std::min(t.maxCoeff(), 0.99f);
            if (sampler->next1D() > p || p <= 0.f) {
                break;
            }
            t /= p;

            /* sample brdf */
            BSDFQueryRecord bRec(its.toLocal(-currRay.d));
            bRec.p = its.p;
            bRec.uv = its.uv;
            Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= brdf;

            /* recursion */
            currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
        }
        
        /* multiple photon passes */
        for (int i = 0; i < m_iterCount; i ++) {
            /* progressive radicance estimate */
            std::vector<uint32_t> results;
            m_photonMaps[i]->search(hp.x, // lookup position
		            m_photonRadius,   // search radius
		            results);
            if (results.size() == 0)
                continue;
            // float density = (hp.N + results.size()) * INV_PI / (hp.R * hp.R);

            /* post-process: radius reduction */
            int N = (int) (hp.N + m_alpha * results.size());
            float dR = hp.R - hp.R * sqrt((float) N / (float) (hp.N + results.size()));
            float R = hp.R - dR;

            /* post-process: flux correction */
            Color3f taoM = 0;
            for (uint32_t j : results) {
                const Photon &photon = (*m_photonMaps[i])[j];
                BSDFQueryRecord bRec(hp.w.normalized(), Frame(hp.n).toLocal(photon.getDirection()), ESolidAngle);
                bRec.p = hp.x;
                bRec.uv = hp.uv;
                taoM += hp.t * hp.bsdf->eval(bRec) * photon.getPower();
            }

            hp.tao = (hp.tao + taoM) * N / (hp.N + results.size());
            hp.N = N;
            hp.R = R;
        }

        if (hp.N == 0 || hp.R == 0) {
            return color;
        } else {
            return (color + INV_PI * hp.tao / (m_emittedPhotonCount * hp.R * hp.R));
        }
    }

    virtual std::string toString() const override {
        return tfm::format(
            "PhotonMapper[\n"
            "  photonCount = %i,\n"
            "  photonRadius = %f\n"
            "]",
            m_photonCount,
            m_photonRadius
        );
    }
private:
    /* 
     * Important: m_photonCount is the total number of photons deposited in the photon map,
     * NOT the number of emitted photons. You will need to keep track of those yourself.
     */ 
    int m_photonCount;
    float m_photonRadius;
    int m_iterCount;
    int m_emittedPhotonCount;
    std::vector<int> m_emittedPhotonCounts;
    std::vector<std::unique_ptr<PhotonMap>> m_photonMaps;
    float m_alpha;
};

NORI_REGISTER_CLASS(ProgressivePhotonMapper, "progressive_photonmapper");
NORI_NAMESPACE_END
