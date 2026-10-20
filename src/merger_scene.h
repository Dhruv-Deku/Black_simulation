#pragma once

#include "audio.h"
#include "camera.h"
#include "gl_util.h"
#include "merger.h"
#include "merger_sim.h"
#include "settings.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

class Renderer;
struct ImDrawList;
struct ImFont;

// The neutron star merger as a scene: the simulation (src/merger_sim.*), the clock that
// carries it from the last orbits (milliseconds) to the kilonova (days), a camera that
// follows the action, the two ways of looking at it, and its sound.
//
//   Cinematic   the scene ray traced as a film would show it (shaders/merger_view.frag), the
//               camera circling it: the stars lensing each other and the stars behind, the
//               spiral of gravitational waves, the collision's flash and the cloud of debris
//               it throws out, the merged remnant, the black hole, the jet and its cocoon,
//               the kilonova's blue and red ejecta;
//   Simulation  the simulation's matter seen from above as numerical relativity papers show
//               it (shaders/merger_map.frag), with the gravitational wave and the clock;
//   Side by side, both at once;
//   Film        the reference short film's shots, replayed with its timing, framing and
//               colours (src/film.h, shaders/merger_film.frag), captions and all.
//
// Sound: the gravitational wave played at its true pitch (it is in the audible band: a
// few hundred Hz rising to kilohertz), following the picture's slow motion.
class MergerScene {
public:
    enum View { Cinematic = 0, Simulation = 1, SideBySide = 2, Film = 3 };

    bool init(const std::filesystem::path& shaderDir, std::string& error);
    void shutdown();
    void hotReload();
    std::string errors() const;

    // Builds the stars and points the camera at them.
    void start(Camera& camera);
    void stop();
    bool active() const { return m_active; }

    // Advances by wallDt seconds of viewing (at the tempo of src/merger.h), with at most
    // stepBudget hydrodynamic steps. live: play the sound now (not while recording).
    void update(double wallDt, bool paused, Camera& camera, int stepBudget, bool live);
    // Draws the frame into the renderer's output image, overlays included.
    void render(Renderer& renderer, const Settings& settings, const Camera& camera, int width, int height,
                float renderScale);
    // The panel's controls: 1 if (re)starting is asked for, 2 to leave the merger.
    int drawControls(Camera& camera);

    // The sound made during the last update (mono, kSampleRate), for videos.
    const std::vector<int16_t>& lastAudio() const { return m_lastAudio; }
    static constexpr int kSampleRate = 44100;

    // The film's clock: seconds into the reference film (it loops).
    double filmTime() const { return m_filmTime; }
    void setFilmTime(double t);
    // The bold italic face the film's captions are set in (Segoe UI Bold Italic).
    void setCaptionFont(ImFont* font) { m_captionFont = font; }

    double sinceMerger() const;  // physical seconds (negative before)
    std::string phaseName() const;
    std::string describe() const;  // a line for logs and screenshots
    const MergerSim& sim() const { return m_sim; }

    // Options (the panel and the command line set them).
    MergerSim::Params params;
    int view = Cinematic;  // the ray-traced 3D scene; the film view only replays the reference
    bool captions = true;
    bool waves = true;
    bool kilonovaOn = true;
    bool cloudOn = true;   // the debris cloud the collision throws out
    bool jetOn = true;
    bool sound = true;
    bool autoCamera = true;
    bool filmColor = true;
    float speed = 1.0f;     // times the tempo
    float pitch = 1.0f;     // sound: times the wave's true frequency
    float volume = 0.6f;
    merger::Jet jet;
    merger::Tempo tempo;
    merger::DebrisCloud debris;

private:
    struct Overlay;
    void directCamera(double wallDt, Camera& camera);
    void synthesize(double wallDt, bool live, bool paused);
    void buildWaveTexture(double window);
    void drawCinematic(Renderer& renderer, const Settings& settings, const Camera& camera, int x, int y, int w, int h);
    void drawMap(int x, int y, int w, int h);
    void drawFilm(Renderer& renderer, int w, int h);
    void fillFilmCaptions(ImDrawList& list, int width, int height);
    void synthesizeFilm(double wallDt, bool live, bool paused);
    void drawOverlays(int width, int height, int splitX);
    void fillOverlays(ImDrawList& list, int width, int height, int splitX);
    double kilonovaRamp(double t) const;
    // The debris cloud as drawn, t after the merger: its optical depth through its radius
    // (0 before the merger and once the kilonova has taken over), and the collision's flash.
    double cloudDepth(double t) const;
    double mergeFlash(double t) const;
    double burstTime() const;
    double ejectaFront(double t) const;  // the fastest ejecta's edge t after the merger (units)
    // The clock stops here: the kilonova has faded a hundredfold and turned infrared.
    static constexpr double kEndTime = 12.0 * 86400.0;
    double scale(const Camera& camera) const;  // simulation units per render unit

    bool m_active = false;
    MergerSim m_sim;
    std::vector<merger::Kilonova> m_kilonova;
    ShaderProgram m_view, m_map, m_film, m_noiseBake;
    GLuint m_waveTexture = 0, m_skyTexture = 0;
    // Tileable noise for the clouds' billows (shaders/merger_noise.comp), baked at the start.
    GLuint m_noiseTexture = 0;
    static constexpr int kNoiseDims = 128;
    void bakeNoise();
    static constexpr int kWaveSamples = 1024;
    double m_wavePeak = 1e-30;
    double m_waveWindow = 1.0;
    bool m_skyUploaded = false;

    // The camera's director.
    double m_lastR = -1.0, m_lastTheta = -1.0, m_lastPhi = -1.0;
    bool m_autoAngles = true;
    double m_viewHalf = 20.0;   // the simulation view's half height (units)
    double m_boxHalf = 20.0;    // the drawing grid's half size
    double m_wallTime = 0.0;
    double m_orbitTime = 0.0;   // seconds the camera has been circling by itself
    // The line through the stars' centres as they touch (angle about the orbit's axis): the
    // tidal tails leave from its ends.
    double m_contactPhase = 0.0;

    // Sound.
    AudioStream m_audio;
    std::vector<int16_t> m_lastAudio;
    double m_audioPhase = 0.0, m_audioFreq = 0.0, m_audioAmp = 0.0;
    double m_sampleCarry = 0.0;
    double m_boomAge = 1e9;     // seconds of viewing since the gamma-ray burst broke out
    bool m_burstSeen = false;
    uint32_t m_noise = 12345;
    double m_lowpass = 0.0;

    // The film.
    double m_filmTime = 0.0;
    double m_filmAudioTime = 0.0;   // where the film's sound has got to
    double m_filmChirpPhase = 0.0;
    double m_boomLow = 0.0, m_rumbleLow = 0.0, m_whooshLow = 0.0;
    ImFont* m_captionFont = nullptr;

    // The strain trace on the simulation view: recent h+ (face-on), newest last.
    std::vector<float> m_trace;
};
