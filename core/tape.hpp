#ifndef TAPE_HPP
#define TAPE_HPP

namespace quadra
{
    class TapeContext
    {
    public:
        TapeContext() { start_recording(); }
        ~TapeContext() { stop_recording(); }
    };
} // namespace pelagia
#endif // TAPE_HPP