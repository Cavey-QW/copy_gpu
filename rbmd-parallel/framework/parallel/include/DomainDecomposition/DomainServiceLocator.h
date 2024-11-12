#pragma once

#include <memory>
#include "Domain.h"

class DomainServiceLocator {
private:
    DomainServiceLocator() {

    }

    ~DomainServiceLocator() = default;

    DomainServiceLocator(const DomainServiceLocator &);

    DomainServiceLocator &operator=(const DomainServiceLocator &);
    std::shared_ptr<Domain> _domain = nullptr;
public:
    static DomainServiceLocator &Instance();

    std::shared_ptr<Domain> GetDomain();

    void SetDomain(Domain*  domain);

};

