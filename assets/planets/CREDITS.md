# Planet maps

Public-domain imagery from NASA, used to draw the Earth and Jupiter when they sit at the
center (src/renderer.cpp, `loadPlanetMaps`; shaders/central_surface.glsl, "Real maps").

| File | Source |
|---|---|
| `earth_day.jpg` | Blue Marble Next Generation, December 2004, with topography and bathymetry (5400 x 2700). NASA Earth Observatory, Reto Stockli. https://visibleearth.nasa.gov/images/73909 |
| `earth_clouds.jpg` | Blue Marble cloud cover, combined (2048 x 1024). NASA Earth Observatory. https://visibleearth.nasa.gov/images/57747 |
| `earth_night.jpg` | Black Marble 2016, city lights (3600 x 1800). NASA Earth Observatory, Joshua Stevens, Miguel Roman. https://earthobservatory.nasa.gov/features/NightLights |
| `jupiter.jpg` | Cassini global map of Jupiter, PIA07782 (3601 x 1801). NASA/JPL/Space Science Institute. https://photojournal.jpl.nasa.gov/catalog/PIA07782 |

All are equirectangular, east to the right, north up. NASA imagery is not copyrighted
(see https://www.nasa.gov/nasa-brand-center/images-and-media/); credit is given above.
