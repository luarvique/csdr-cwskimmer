#include "bufmodule.hpp"

using namespace Csdr;

template <typename T, typename U>
BufferedModule<T, U>::BufferedModule(Csdr::Module<T, U> *module, size_t bufferSize)
{
  this->module = module;
  this->buffer = new Csdr::Ringbuffer<T>(bufferSize);
  this->reader = new Csdr::RingbufferReader<T>(this->buffer);
  this->module->setReader(this->reader);
}

template <typename T, typename U>
BufferedModule<T, U>::~BufferedModule()
{
  delete module;
  delete reader;
  delete buffer;
}

template <typename T, typename U>
void BufferedModule<T, U>::connect(BufferedModule<U, U> *output)
{
  setWriter(output->buf());
}

template <typename T, typename U>
void BufferedModule<T, U>::processAll()
{
  while(canProcess()) process();
}

namespace Csdr
{
  template class BufferedModule<float, unsigned char>;
  template class BufferedModule<unsigned char, unsigned char>;
}
