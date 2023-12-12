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
    Color3f wgt;    // pixel weight
    float R;        // current photon radius
    int N;          // accumulated photon count
    Color3f tao;    // accumulated reflected flux
};

class ProgressivePhotonMapper : public Integrator {
public:
    /// Photon map data structure
    typedef PointKDTree<Photon> PhotonMap;
    int m_emittedPhotonCount = 0;

    ProgressivePhotonMapper(const PropertyList &props) {
        /* Lookup parameters */
        m_photonCount  = props.getInteger("photonCount", 1000000);
        m_photonRadius = props.getFloat("photonRadius", 0.0f /* Default: automatic */);
        m_iterCount = props.getInteger("interCount", 10);
        m_alpha = props.getFloat("alpha", 0.7);
    }

    virtual void preprocess(const Scene *scene) override {
        cout << "Gathering " << m_photonCount << " photons .. ";
        cout.flush();

        /* Create a sample generator for the preprocess step */
        Sampler *sampler = static_cast<Sampler *>(
            NoriObjectFactory::createInstance("independent", PropertyList()));

        /* Allocate memory for the photon map */
        m_photonMap = std::unique_ptr<PhotonMap>(new PhotonMap());
        m_photonMap->reserve(m_photonCount);

		/* Estimate a default photon radius */
		if (m_photonRadius == 0)
			m_photonRadius = scene->getBoundingBox().getExtents().norm() / 500.0f;
	

		/* How to add a photon?
		 * m_photonMap->push_back(Photon(
		 *	Point3f(0, 0, 0),  // Position
		 *	Vector3f(0, 0, 1), // Direction
		 *	Color3f(1, 2, 3)   // Power
		 * ));
		 */

        /* ray tracing pass*/

        /* sample and trace ray for each pixel */
        const Camera *camera = scene->getCamera();
        const Vector2i size = camera->getOutputSize();

        for (int x = 0; x < size.x(); x ++) {
            for (int y = 0; y < size.y(); y ++) {
                Point2f pixelSample = Point2f(x, y) + sampler->next2D();
                Point2f apertureSample = sampler->next2D();

                Ray3f ray;
                Color3f wgt = camera->sampleRay(ray, pixelSample, apertureSample);

                /* store all hit points where the surface has a non-specular component in the BRDF */
                Intersection its;
                Ray3f currRay = ray;
                Color3f t = 1.f;
                Color3f color = 0;
                while (true) {
                    if (!scene->rayIntersect(currRay, its)) {
                        break;
                    }

                    /* add hit points for diffuse surfaces */
                    if (its.mesh->getBSDF()->isDiffuse()) {
                        HitPoint hp;
                        hp.x = its.p;
                        hp.n = its.shFrame.n;
                        hp.w = currRay.d;
                        hp.bsdf = its.mesh->getBSDF();
                        hp.xy = Point2i(x, y);
                        hp.wgt = wgt;
                        hp.R = m_photonRadius;
                        hp.N = 0;
                        hp.tao = 0;

                        m_hitpoints.push_back(hp);
                        break;
                    }

                    /* Russian roulette */
                    float p = std::min(t.maxCoeff(), 0.99f);
                    if (sampler->next1D() > p || p <= 0.f) {
                        break;
                    }
                    t /= p;

                    /* sample brdf */
                    BSDFQueryRecord bRec(its.shFrame.toLocal(-currRay.d));
                    bRec.p = its.p;
                    bRec.uv = its.uv;
                    Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                    t *= brdf;

                    /* recursion */
                    currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
                }
            }
        }
    }

    virtual Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
    	
		/* How to find photons?
		 * std::vector<uint32_t> results;
		 * m_photonMap->search(Point3f(0, 0, 0), // lookup position
		 *                     m_photonRadius,   // search radius
		 *                     results);
		 *
		 * for (uint32_t i : results) {
		 *    const Photon &photon = (*m_photonMap)[i];
		 *    cout << "Found photon!" << endl;
		 *    cout << " Position  : " << photon.getPosition().toString() << endl;
		 *    cout << " Power     : " << photon.getPower().toString() << endl;
		 *    cout << " Direction : " << photon.getDirection().toString() << endl;
		 * }
		 */

        Intersection its;
        Ray3f currRay = _ray;
        Color3f t = 1.f;
        Color3f color = 0;
        HitPoint hp;

        /* store hit point */ 
        while (true) {
            if (!scene->rayIntersect(currRay, its)) {
                break;
            }

            /* add hit points for diffuse surfaces */
            if (its.mesh->getBSDF()->isDiffuse()) {
                hp.x = its.p;
                hp.n = its.shFrame.n;
                hp.w = currRay.d;
                hp.bsdf = its.mesh->getBSDF();
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
            BSDFQueryRecord bRec(its.shFrame.toLocal(-currRay.d));
            bRec.p = its.p;
            bRec.uv = its.uv;
            Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            t *= brdf;

            /* recursion */
            currRay = Ray3f(bRec.p, its.toWorld(bRec.wo));
        }
        
        /* multiple photon passes */
        int m_emittedPhotonCount = 0;
        for (int i = 0; i < m_iterCount; i ++) {
            /* trace photon */
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
                        m_photonMap->push_back(Photon(
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

            /* Build the photon map */
            m_photonMap->build();

            /* progressive radicance estimate */
            std::vector<uint32_t> results;
            m_photonMap->search(hp.x, hp.R, results);
            float density = (hp.N + results.size()) * INV_PI / hp.R;

            /* post-process: radius reduction */
            int nextN = (int) hp.N + m_alpha * results.size();
            float dR = hp.R - hp.R * sqrt(nextN / (hp.N + results.size()));

            /* post-process: flux correction */
            Color3f nextFlux = 0;
            for (uint32_t i : results) {
                const Photon &photon = (*m_photonMap)[i];
                nextFlux += photon.getPower();
            }
            nextFlux = (hp.tao + nextFlux) * nextN / (hp.N + results.size());

            hp.N = nextN;
            hp.R -= dR;
            hp.tao = nextFlux;

            m_photonMap->clear();
        }

        /* radiance evaluation */
		return INV_PI * hp.tao / (m_emittedPhotonCount * hp.R * hp.R);
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
    std::unique_ptr<PhotonMap> m_photonMap;
    std::vector<HitPoint> m_hitpoints;
    float m_alpha;
};

NORI_REGISTER_CLASS(ProgressivePhotonMapper, "progressive_photonmapper");
NORI_NAMESPACE_END
