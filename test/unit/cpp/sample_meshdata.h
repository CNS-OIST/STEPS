#pragma once

double COORDS[][3] = {
    {0.0, 0.0, 0.0},
    {1.0, 0.0, 0.0},
    {0.5, 0.866025403784, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 0.0, 1.0},
    {0.5, 0.866025403784, 1.0},
};
steps::vertex_id_t::value_type TETINDICES[][4] = {
    {0uL, 1uL, 2uL, 3uL},
    {1uL, 4uL, 3uL, 2uL},
    {2uL, 4uL, 5uL, 3uL},
};
