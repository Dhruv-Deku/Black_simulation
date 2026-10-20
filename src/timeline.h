#pragma once

// Scrubbing back and forth through what has happened.
//
// The simulation is deterministic: from the same state, the same frame steps do exactly
// the same thing, bit for bit (as saving and loading scenes shows; src/scene_file.h).
// So the timeline keeps checkpoints, full snapshots of the scene taken every so often,
// and the step every frame took between them. To show an earlier moment it restores the
// last checkpoint before it and replays the logged steps up to it, without drawing them.
// Played on, the replay retraces exactly what happened, up to the present.
//
// Changing a setting that changes what happens (gravity between bodies, the gas physics,
// the central object) is logged with the next step and replayed with it. What the user
// does to the scene itself (adding a body, a starquake) is not: a checkpoint is taken
// right after it, and kept, so that no replay runs across it from before. Either one,
// done while looking at the past, starts a new future: the old one is dropped.
//
// Checkpoints are thinned out as memory fills (every other one goes, doubling the
// spacing), so the whole history stays reachable, at the cost of longer replays.

#include "scene_file.h"

#include <cstddef>
#include <string>
#include <vector>

class Timeline {
public:
    // One frame of simulation.
    struct Step {
        double dt = 0.0;       // simulated time (M)
        float realDt = 0.0f;   // the frame's real time, in seconds (the camera's easing uses it)
        std::string settings;  // if not empty: the simulation settings from this step on (scene::simulationSettings)
    };
    struct Checkpoint {
        size_t step = 0;       // the steps before it lead up to it
        bool pinned = false;   // taken right after something the user did: never thinned out
        std::string scene;     // the scene as JSON (src/scene_file.h)
        scene::Blob blob;      // ... and its binary parts
        int gpu = -1;          // the app's copy of the gas on the GPU (-1: none)
        size_t gpuBytes = 0;
    };

    void clear();
    bool empty() const { return m_checkpoints.empty(); }
    // Times of the first checkpoint, of the latest step, and of the scene now.
    double start() const { return m_times.empty() ? 0.0 : m_times.front(); }
    double end() const { return m_times.empty() ? 0.0 : m_times.back(); }
    double time() const { return m_times.empty() ? 0.0 : m_times[m_position]; }
    // Steps are counted from the first checkpoint.
    size_t position() const { return m_position; }
    size_t length() const { return m_steps.size(); }
    bool atEnd() const { return m_position == m_steps.size(); }

    // --- Recording (at the end) ------------------------------------------------------
    // A checkpoint of the scene as it is now; the first one starts the timeline at `time`.
    // One already taken at this step is replaced: its GPU copy is returned to free (or -1).
    int addCheckpoint(Checkpoint checkpoint, double time);
    // A step just taken, which brought the scene to timeAfter.
    void addStep(Step step, double timeAfter);
    // Time of the latest checkpoint.
    double lastCheckpointTime() const;
    // The simulation settings in force now (as last logged or checkpointed).
    const std::string& settings() const { return m_settings; }
    void setSettings(const std::string& settings) { m_settings = settings; }

    // --- Seeking ----------------------------------------------------------------------
    // The latest checkpoint at or before step `target`.
    const Checkpoint& checkpointFor(size_t target) const;
    const Step& step(size_t index) const { return m_steps[index]; }
    double timeAt(size_t position) const { return m_times[position]; }
    // The last position whose time is at or before `time` (the first if none is).
    size_t positionAt(double time) const;
    // The scene is now at this position (after restoring and replaying up to it).
    void setPosition(size_t position);

    // A new future from here: drops every step and checkpoint after the current position.
    // Their GPU copies are returned for the app to free.
    std::vector<int> truncate();

    // --- Memory -------------------------------------------------------------------------
    size_t bytes() const;     // checkpoints and log on the CPU
    size_t gpuBytes() const;  // the gas copies
    size_t checkpointCount() const { return m_checkpoints.size(); }
    size_t pinnedCount() const;
    // Drops every other unpinned checkpoint (keeping the first and the latest) while over
    // either budget; returns the GPU copies to free.
    std::vector<int> thin(size_t cpuBudget, size_t gpuBudget);

private:
    std::vector<Step> m_steps;
    std::vector<double> m_times;  // m_times[i]: the time after i steps (m_times[0]: the first checkpoint's)
    std::vector<Checkpoint> m_checkpoints;  // in order of step
    size_t m_position = 0;
    std::string m_settings;
};
