#pragma once

// The telescope's view: what the Event Horizon Telescope would see of the black hole
// (src/eht.h), the model image next to the same image blurred to the array's resolution,
// at the real angular scale, with the shadow's edge, the measured ring and the flux.
// Images are made on worker threads (a quick preview first) whenever an option changes.

#include "eht.h"

#include <glad/gl.h>

#include <atomic>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

class EhtView {
public:
    enum Target { M87, SgrA, ThisHole };
    struct Options {
        int target = M87;
        double spin = 0.5;
        double inclinationDeg = 163.0;
        double frequencyGHz = 230.0;
        double temperature = 1.5e12;  // T0, K at r = M
        double beta = 10.0;
        double height = 0.5;        // the flow's H/R
        int size = 160;             // pixels across
        double beamUas = 20.0;      // the array's resolution (FWHM)
        double massSolar = 6.5e9;   // for ThisHole
        bool curve = true;          // draw the shadow's edge (critical curve)
        bool operator==(const Options&) const = default;
    };

    ~EhtView();
    // Starts a new image when the options change; takes in finished ones (call each frame,
    // on the thread with the GL context).
    void update(const Options& options);
    // Draws the view over the window right of x = left (pixels), above the bottom margin.
    void draw(float left, float bottomMargin);
    bool busy() const { return m_running.load(); }
    // Waits for the image being made (for screenshots).
    void finish();
    void shutdown();

    // The source's numbers, for any target.
    static const char* targetName(int target);

private:
    struct Result {
        Options options;
        eht::Image model, blurred;
        std::vector<eht::Point> curve;
        double thetaG = 0.0;   // uas per M
        double n0 = 0.0;
        double flux = 0.0;     // Jy
        double ringUas = 0.0;  // diameter of the blurred ring
        double shadowUas = 0.0;  // mean diameter of the critical curve
        double brightest = 0.0;  // brightness temperature of the brightest pixel, K (model)
        bool preview = false;
        double seconds = 0.0;
    };
    void work(Options options);
    void upload(const Result& r);
    void stopWorker();

    Options m_options;
    bool m_started = false;
    std::thread m_worker;
    std::atomic<bool> m_cancel{false}, m_running{false};
    std::mutex m_mutex;
    Result m_pending;
    bool m_hasPending = false;
    Result m_shown;
    bool m_hasShown = false;
    GLuint m_textures[2] = {0, 0};
    std::string m_status;
};
