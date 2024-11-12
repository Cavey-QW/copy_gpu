#include <stdexcept>
#include "DomainServiceLocator.h"

DomainServiceLocator &DomainServiceLocator::Instance() {
    static DomainServiceLocator instance;
    return instance;
}

std::shared_ptr<Domain> DomainServiceLocator::GetDomain() {
    if (_domain == nullptr)
        throw std::runtime_error("Domain instance is not set");
    return _domain;
}

void DomainServiceLocator::SetDomain(Domain* domain) {
    _domain.reset(domain);
}
