#ifndef BUFMODULE_HPP
#define BUFMODULE_HPP

#include "csdr/module.hpp"
#include "csdr/ringbuffer.hpp"

namespace Csdr {

template <typename T, typename U> class BufferedModule: public Module<T, U>
{
  public:
    BufferedModule(Module<T, U> *module, size_t bufferSize);
    ~BufferedModule();

    bool canProcess() override { return(module->canProcess()); }
    void process() override { module->process(); }

    void wait(std::unique_lock<std::mutex>& lock) override { module->wait(lock); }
    void unblock() override { module->unblock(); }
    void setWriter(Writer<U>* writer) override { module->setWriter(writer); }
    void setReader(Reader<T>* reader) override { /* do nothing */ }

    void connect(BufferedModule<U, U> *output);
    void processAll();

    Module<T, U>  *mod() { return(module); }
    Ringbuffer<T> *buf() { return(buffer); }
    RingbufferReader<T> *rdr() { return(reader); }

  private:
    Module<T, U> *module;
    Ringbuffer<T> *buffer;
    RingbufferReader<T> *reader;
};

}

#endif // BUFMODULE_HPP
