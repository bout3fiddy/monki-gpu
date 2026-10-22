"""Common lifecycle for interchangeable transport implementations."""

from abc import ABC, abstractmethod
from types import TracebackType
from typing import Self

from ..output import SimulationResult


class Backend(ABC):
    """Strategy interface consumed by the shared simulation driver.

    Implementations own prepared tables and reusable execution memory. The driver
    owns rounds, convergence, and normalization so statistics are shared across devices.

    CPU and Metal accept ReferenceCase and BatchControls, and return
    RawTallies. Prepared metadata supplies labels
    and optical fields; native physics routines receive explicit named arrays.
    """

    metadata: SimulationResult | None = None
    provenance: dict = {}

    @abstractmethod
    def prepare(self, scene: object) -> None:
        """Validate and retain scene tables before running batches.

        Repreparing replaces prior scene state. Failure must leave resources safe
        to close; preparation must succeed before run_batch can execute.
        """

        raise NotImplementedError

    @abstractmethod
    def run_batch(self, controls: object) -> object:
        """Return owned raw tallies for the requested photon range, without normalization.

        Returned data must remain valid after the next batch or close operation.
        The simulation driver owns convergence decisions.
        """

        raise NotImplementedError

    def run_batches(self, controls: list) -> list:
        """Return one owned result per requested round, preserving order and statistics.

        Backends may schedule independent rounds together. The default retains the
        sequential interface, so CPU and future backends need no special scheduler.
        """

        return [self.run_batch(request) for request in controls]

    @abstractmethod
    def close(self) -> None:
        """Release resources; safe before preparation and on repeated calls."""

        raise NotImplementedError

    def __enter__(self) -> Self:
        """Return this strategy without implicitly allocating hardware resources."""

        return self

    def __exit__(
        self,
        exception_type: type[BaseException] | None,
        exception: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        """Release resources even if preparation or simulation raises an error."""

        self.close()
