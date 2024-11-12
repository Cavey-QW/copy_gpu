#pragma once

#include <vector>
#include <cstddef>
#include <mpi.h>
#include "LogTool.h"
#include "Atom.h"
// do not uncomment the if, it will break halo copies of the kddecomposition!
//#if (not defined(NDEBUG))
#define LS1_SEND_UNIQUE_ID_FOR_HALO_COPIES
//#pragma message "Compilation info: Unique IDs of Halo-Molecules are always present."
//#endif

/**
    这个类使得在发送HALO分子时只发送位置数据
    而在发送LEAVING分子时发送所有数据
    它首先发送两个整数,告知要发送多少个LEAVING分子和多少个HALO分子
    然后发送所有数据
    内部将所有数据转换为CHAR     ----- 字节 char = byte 指针其实是一段空间
    TODO: 测试在大端/小端架构上由于CHAR转换会如何工作           ------------ CommunicationBuffer 和 AppendInt怎么说
    存储两个无符号长整型,然后是LEAVING分子,最后是HALO分子
    */
class CommunicationBuffer {
public:
	CommunicationBuffer() {
		clear();
	}
	size_t getDynamicSize();

	void clear();

	void resizeForAppendingLeavingMolecules(unsigned long numMols);
	void resizeForAppendingHaloMolecules(unsigned long numMols);
        void resizeForAppendingForceMolecules(unsigned long numMols);

	unsigned char * getDataForSending();
	size_t getNumElementsForSending();
	void resizeForRawBytes(unsigned long numBytes);

	// write
	void addLeavingMolecule(size_t indexOfMolecule, const Atom& m);
	void addHaloMolecule(size_t indexOfMolecule, const Atom& m);
    void addForceMolecule(size_t indexOfMolecule, const Atom& m);

	// read
	void readLeavingMolecule(size_t indexOfMolecule, Atom& m) const;
	void readHaloMolecule(size_t indexOfMolecule, Atom& m) const;
	void readForceMolecule(size_t indexOfMolecule, Atom& m) const;

	void resizeForReceivingMolecules(unsigned long& numLeaving, unsigned long& numHalo);
	void resizeForReceivingMolecules(unsigned long& numForces);

	size_t getNumHalo() const {
		return _numHalo;
	}

	size_t getNumLeaving() const {
		return _numLeaving;
	}

        size_t getNumForces() const {
            return _numForces;
        }

	static MPI_Datatype getMPIDataType() {
		return MPI_CHAR;
	}

private:
	static size_t _numBytesHalo;
	static size_t _numBytesLeaving;
        static size_t _numBytesForces; // where is this set?

	enum class ParticleType_t {HALO=0, LEAVING=1, FORCE=3};
	size_t getStartPosition(ParticleType_t type, size_t indexOfMolecule) const;

	/**
	 * @return the next index for writing
	 */
	template<typename T>
	size_t emplaceValue(size_t indexInBytes, T passByValue);

	template<typename T>
	size_t readValue(size_t indexInBytes, T& passByReference) const;

	typedef unsigned char byte_t;
	std::vector<byte_t> _buffer;
	size_t _numLeaving, _numHalo, _numForces;
};

template<typename T>
inline size_t CommunicationBuffer::emplaceValue(size_t indexInBytes, T passByValue) {
	const size_t numBytesOfT = sizeof(T);
	size_t ret = indexInBytes + numBytesOfT;

	//assert(_buffer.size() >= ret);
    if (_buffer.size() < ret){
        GlobalLogger::error("_buffer.size() < ret in CommunicationBuffer::emplaceValue");
    }


	const byte_t * pointer = reinterpret_cast<byte_t *> (&passByValue);
	for (size_t i = 0; i < numBytesOfT; ++i) {
		_buffer[indexInBytes + i] = pointer[i];
	}

	return ret;
}

template<typename T>
inline size_t CommunicationBuffer::readValue(size_t indexInBytes, T& passByReference) const {
	const size_t numBytesOfT = sizeof(T);
	size_t ret = indexInBytes + numBytesOfT;

    // assert(_buffer.size() >= ret);
    if (_buffer.size() < ret){
        GlobalLogger::error("_buffer.size() < ret in CommunicationBuffer::readValue");
    }

	byte_t * pointer = reinterpret_cast<byte_t *> (&passByReference);
	for (size_t i = 0; i < numBytesOfT; ++i) {
		pointer[i] = _buffer[indexInBytes + i];
	}

	return ret;
}
