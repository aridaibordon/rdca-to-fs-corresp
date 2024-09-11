from dataclasses import dataclass


@dataclass
class EnergyLevel:
    n: int
    l: int
    s: str

    @property
    def spec_notation(self) -> str:
        """Spectroscopic notation"""
        return f"{self.n}{orb_to_sn(self.l)}{self.s}"

    @property
    def deg(self) -> int:
        """Energy level degeneracy"""
        if self.s == "-":
            return 2 * self.l
        return 2 * (self.l + 1)

    def get_subshell(self) -> str:
        return f"{self.n}{orb_to_sn(self.l)}"

    def __repr__(self) -> str:
        return self.spec_notation


@dataclass
class RDCAConfig:
    config: list[tuple[EnergyLevel, int]]

    def __post_init__(self):
        if not self.is_valid():
            raise ValueError(f"{self.config} is not a valid configuration")

    @property
    def ncomplex(self) -> str:
        data = {}
        for lvl, ocu in self.config:
            subshell_label = lvl.n
            if subshell_label not in data.keys():
                data[subshell_label] = ocu
                continue
            data[subshell_label] += ocu

        return ".".join([f"{key}*{data[key]}" for key in data.keys()])

    @property
    def sname(self) -> str:
        data = {}
        for lvl, ocu in self.config:
            subshell_label = lvl.get_subshell()
            if subshell_label not in data.keys():
                data[subshell_label] = ocu
                continue
            data[subshell_label] += ocu

        for key in data.copy().keys():
            orb_ind, orb_ocu = key[-1], data[key]

            max_ocu = 4 * orb_to_sn(orb_ind) + 2
            if orb_ocu == max_ocu:
                data.pop(key, None)

        if not data:
            return f"{subshell_label}{ocu}"

        return ".".join([f"{key}{data[key]}" for key in data.keys()])

    @property
    def name(self) -> str:
        name = []
        for lvl, ocu in self.config:
            if lvl.deg == ocu:
                continue
            name.append(f"{lvl.spec_notation}{ocu}")

        if not name:
            return f"{lvl}{ocu}"
        return ".".join(name)

    @property
    def full(self) -> dict:
        return {
            "ncomplex": self.ncomplex,
            "sname": self.sname,
            "name": self.name,
        }

    def is_valid(self) -> str:
        for lvl, ocu in self.config:
            if ocu > lvl.deg:
                return False

        return True

    def __getitem__(self, attr):
        return getattr(self, attr)

    def __repr__(self) -> str:
        return str(self.config)


def orb_to_sn(l: int | str):
    """Convert from quantum orbital number to spectroscopic notation."""
    orb_list = "spdfghiklmnoqrtuvwxy"
    if type(l) == str:
        return orb_list.index(l)
    return orb_list[l]


def abako_to_atomic_config(abako_not: str) -> RDCAConfig:
    """Convert ABAKO notation to an atomic configuration"""
    nmax = 10

    energy_lvls = []
    for n in range(1, nmax + 1):
        for l in range(0, n):
            if l != 0:
                energy_lvls.append(EnergyLevel(n, l, "-"))
            energy_lvls.append(EnergyLevel(n, l, "+"))

    config = [
        (energy_lvls[ind], int(ocu)) for ind, ocu in enumerate(abako_not) if int(ocu)
    ]

    return RDCAConfig(config)
