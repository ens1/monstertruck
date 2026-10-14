//! Directed-rounding interval arithmetic for continuation admission.

#[derive(Clone, Copy, Debug)]
pub(super) struct Interval {
    pub(super) low: f64,
    pub(super) high: f64,
}

impl Interval {
    pub(super) const ZERO: Self = Self {
        low: 0.0,
        high: 0.0,
    };
    pub(super) const ONE: Self = Self {
        low: 1.0,
        high: 1.0,
    };
    pub(super) const fn point(value: f64) -> Self {
        Self {
            low: value,
            high: value,
        }
    }
    pub(super) const fn new(low: f64, high: f64) -> Self { Self { low, high } }
    pub(super) fn valid(self) -> bool {
        self.low.is_finite() && self.high.is_finite() && self.low <= self.high
    }
    fn invalid() -> Self { Self::point(f64::NAN) }
    pub(super) fn add(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::invalid();
        }
        if self.low == 0.0 && self.high == 0.0 {
            return other;
        }
        if other.low == 0.0 && other.high == 0.0 {
            return self;
        }
        Self::new(
            (self.low + other.low).next_down(),
            (self.high + other.high).next_up(),
        )
    }
    pub(super) fn neg(self) -> Self { Self::new(-self.high, -self.low) }
    pub(super) fn sub(self, other: Self) -> Self { self.add(other.neg()) }
    pub(super) fn mul(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::invalid();
        }
        if (self.low == 0.0 && self.high == 0.0) || (other.low == 0.0 && other.high == 0.0) {
            return Self::ZERO;
        }
        let values = [
            self.low * other.low,
            self.low * other.high,
            self.high * other.low,
            self.high * other.high,
        ];
        if values.iter().any(|value| !value.is_finite()) {
            return Self::invalid();
        }
        Self::new(
            values.into_iter().fold(f64::INFINITY, f64::min).next_down(),
            values
                .into_iter()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        )
    }
    pub(super) fn div(self, other: Self) -> Self {
        if !other.valid() || (other.low <= 0.0 && other.high >= 0.0) {
            return Self::invalid();
        }
        self.mul(Self::new(
            (1.0 / other.high).next_down(),
            (1.0 / other.low).next_up(),
        ))
    }
    pub(super) fn hull(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::invalid();
        }
        Self::new(self.low.min(other.low), self.high.max(other.high))
    }
    pub(super) fn intersect(self, other: Self) -> Self {
        if !self.valid() || !other.valid() {
            return Self::invalid();
        }
        Self::new(self.low.max(other.low), self.high.min(other.high))
    }
    pub(super) fn absolute(self) -> f64 {
        if self.valid() {
            self.low.abs().max(self.high.abs())
        } else {
            f64::INFINITY
        }
    }
    pub(super) fn middle(self) -> f64 { self.low * 0.5 + self.high * 0.5 }
    pub(super) fn interior_contains(self, other: Self) -> bool {
        self.valid() && other.valid() && self.low < other.low && other.high < self.high
    }
}

pub(super) type Vector = [Interval; 3];
pub(super) type Parameters = [Interval; 4];
pub(super) type Matrix = [[Interval; 4]; 4];

pub(super) fn dot(a: &Vector, b: &Vector) -> Interval {
    a.iter()
        .zip(b)
        .fold(Interval::ZERO, |sum, (a, b)| sum.add(a.mul(*b)))
}

pub(super) fn distance_bound(a: &Vector, b: &Vector) -> f64 {
    let square = a.iter().zip(b).fold(Interval::ZERO, |sum, (a, b)| {
        let distance = Interval::point(a.sub(*b).absolute());
        sum.add(distance.mul(distance))
    });
    if square.valid() {
        square.high.sqrt().next_up()
    } else {
        f64::INFINITY
    }
}
