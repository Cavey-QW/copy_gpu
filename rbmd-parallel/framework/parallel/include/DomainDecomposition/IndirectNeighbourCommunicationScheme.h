#pragma once
#include "BaseNeighbourCommunicationScheme.h"
#include "DomainDecompSequential.h"


class IndirectNeighbourCommunicationScheme:public BaseNeighbourCommunicationScheme {
private:
protected:
    void initExchangeMoleculesMPI1D(LinkedCell* moleculeContainer,  MessageType msgType,
            bool removeRecvDuplicates, unsigned short d, DomainDecompBase* domainDecomp);

    void finalizeExchangeMoleculesMPI1D(LinkedCell* moleculeContainer, MessageType msgType,
            bool removeRecvDuplicates, unsigned short d, DomainDecompBase* domainDecomp);
    void exchangeMoleculesMPI1D(LinkedCell* moleculeContainer,  MessageType msgType,
            bool removeRecvDuplicates, unsigned short d, DomainDecompBase* domainDecomp);
    void convert1StageTo3StageNeighbours(const std::vector<CommunicationPartner>& commPartners,
            std::vector<std::vector<CommunicationPartner>>& neighbours, HaloRegion& ownRegion, double cutoffRadius);
public:
    explicit IndirectNeighbourCommunicationScheme(ZonalMethod* zonalMethod) :
            BaseNeighbourCommunicationScheme(3, zonalMethod, false) {
    }
    ~IndirectNeighbourCommunicationScheme() override = default;

    void exchangeMoleculesMPI(LinkedCell* moleculeContainer,  MessageType msgType,
            bool removeRecvDuplicates, DomainDecompBase* domainDecomp, bool doHaloPositionCheck=true) override;

    void initCommunicationPartners(double cutoffRadius,
            DomainDecompBase* domainDecomp,
            LinkedCell* moleculeContainer) override;


    std::vector<int> get3StageNeighbourRanks() override {
        std::vector<int> neighbourRanks;
        for (auto & _fullShellNeighbour : _fullShellNeighbours) {
            if (_fullShellNeighbour.isFaceCommunicator()) {
                neighbourRanks.push_back(_fullShellNeighbour.getRank());
            }
        }
        return neighbourRanks;
    }

    void prepareNonBlockingStageImpl(LinkedCell* moleculeContainer,
            unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
            DomainDecompBase* domainDecomp) override;

    void finishNonBlockingStageImpl(LinkedCell* moleculeContainer,
            unsigned int stageNumber, MessageType msgType, bool removeRecvDuplicates,
            DomainDecompBase* domainDecomp) override;



};

