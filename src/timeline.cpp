#include "timeline.h"

#include <algorithm>
#include <cassert>

void Timeline::clear() {
    m_steps.clear();
    m_times.clear();
    m_checkpoints.clear();
    m_position = 0;
    m_settings.clear();
}

int Timeline::addCheckpoint(Checkpoint checkpoint, double time) {
    assert(atEnd());
    if (m_times.empty()) m_times.push_back(time);
    checkpoint.step = m_position;
    // A second checkpoint at the same step (something else the user did before any time
    // passed) replaces the first.
    if (!m_checkpoints.empty() && m_checkpoints.back().step == m_position) {
        const int released = m_checkpoints.back().gpu;
        checkpoint.pinned = checkpoint.pinned || m_checkpoints.back().pinned;
        m_checkpoints.back() = std::move(checkpoint);
        return released;
    }
    m_checkpoints.push_back(std::move(checkpoint));
    return -1;
}

void Timeline::addStep(Step step, double timeAfter) {
    assert(atEnd() && !m_checkpoints.empty());
    m_steps.push_back(std::move(step));
    m_times.push_back(timeAfter);
    m_position = m_steps.size();
}

double Timeline::lastCheckpointTime() const {
    return m_checkpoints.empty() ? -1e300 : m_times[m_checkpoints.back().step];
}

const Timeline::Checkpoint& Timeline::checkpointFor(size_t target) const {
    // Checkpoints are in order of step: the last one not after the target.
    auto it = std::upper_bound(m_checkpoints.begin(), m_checkpoints.end(), target,
                               [](size_t t, const Checkpoint& c) { return t < c.step; });
    return it == m_checkpoints.begin() ? m_checkpoints.front() : *(it - 1);
}

size_t Timeline::positionAt(double time) const {
    if (m_times.empty()) return 0;
    // Times never decrease along the log (steps are never negative).
    auto it = std::upper_bound(m_times.begin(), m_times.end(), time);
    return it == m_times.begin() ? 0 : static_cast<size_t>(it - m_times.begin()) - 1;
}

void Timeline::setPosition(size_t position) {
    m_position = std::min(position, m_steps.size());
}

std::vector<int> Timeline::truncate() {
    std::vector<int> released;
    while (!m_checkpoints.empty() && m_checkpoints.back().step > m_position) {
        if (m_checkpoints.back().gpu >= 0) released.push_back(m_checkpoints.back().gpu);
        m_checkpoints.pop_back();
    }
    m_steps.resize(m_position);
    m_times.resize(m_position + 1);
    return released;
}

size_t Timeline::bytes() const {
    size_t total = m_steps.size() * sizeof(Step) + m_times.size() * sizeof(double);
    for (const Step& s : m_steps) total += s.settings.size();
    for (const Checkpoint& c : m_checkpoints) {
        total += c.scene.size() + c.blob.bytes();
    }
    return total;
}

size_t Timeline::pinnedCount() const {
    return static_cast<size_t>(std::count_if(m_checkpoints.begin(), m_checkpoints.end(),
                                             [](const Checkpoint& c) { return c.pinned; }));
}

size_t Timeline::gpuBytes() const {
    size_t total = 0;
    for (const Checkpoint& c : m_checkpoints) total += c.gpuBytes;
    return total;
}

std::vector<int> Timeline::thin(size_t cpuBudget, size_t gpuBudget) {
    std::vector<int> released;
    while ((bytes() > cpuBudget || gpuBytes() > gpuBudget) && m_checkpoints.size() > 2) {
        // Every other unpinned checkpoint between the first and the latest goes.
        std::vector<Checkpoint> kept;
        kept.reserve(m_checkpoints.size());
        bool drop = false;
        for (size_t i = 0; i < m_checkpoints.size(); ++i) {
            Checkpoint& c = m_checkpoints[i];
            const bool candidate = i > 0 && i + 1 < m_checkpoints.size() && !c.pinned;
            if (candidate && drop) {
                if (c.gpu >= 0) released.push_back(c.gpu);
            } else {
                kept.push_back(std::move(c));
            }
            if (candidate) drop = !drop;
        }
        if (kept.size() == m_checkpoints.size()) break;  // only pinned ones are left
        m_checkpoints = std::move(kept);
    }
    return released;
}
