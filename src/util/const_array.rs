use core::alloc::Layout;

use crate::util;


pub struct ConstArray<T>
{
    ptr: *mut T,
    len: usize
}

impl<T> ConstArray<T>
{
    pub const fn new() -> Self
    {
        Self {
            ptr: core::ptr::null_mut(),
            len: 0
        }
    }

    pub const fn len(&self) -> usize
    {
        self.len
    }

    const fn layout(len: usize) -> Layout
    {
        let layout = unsafe {
            Layout::for_value_raw(core::ptr::slice_from_raw_parts(core::ptr::null::<T>(), len))
        };
        match layout.align_to(layout.align().next_power_of_two())
        {
            Ok(l) => l,
            Err(_) => panic!("Invalid align.")
        }
    }

    pub const fn destruct(mut self)
    {
        self.clear();
    }

    pub const fn clear(&mut self)
    {
        util::panic_in_rt();
        let old_layout = ConstArray::<T>::layout(self.len);
        if !self.ptr.is_null()
        {
            unsafe {
                core::intrinsics::const_deallocate(self.ptr.cast(), old_layout.size(), old_layout.align());
            }
        }
        self.ptr = core::ptr::null_mut();
        self.len = 0;
    }

    pub const fn push(&mut self, value: T)
    {
        util::panic_in_rt();
        let new_len = self.len + 1;
        let old_layout = Self::layout(self.len);
        let new_layout = Self::layout(new_len);
        self.ptr = unsafe {
            let new_ptr = core::intrinsics::const_allocate(new_layout.size(), new_layout.align()).cast::<T>();
            if !self.ptr.is_null()
            {
                core::ptr::copy_nonoverlapping(self.ptr, new_ptr, self.len);
                core::intrinsics::const_deallocate(self.ptr.cast(), old_layout.size(), old_layout.align());
            }
            let dst = new_ptr.add(self.len);
            dst.write(value);
            new_ptr
        };
        self.len = new_len
    }

    pub const fn as_slice(&self) -> &[T]
    {
        let &Self {ptr, len} = self;
        if ptr.is_null()
        {
            return &[]
        }
        unsafe {
            core::slice::from_raw_parts(ptr, len)
        }
    }

    pub const fn leak(self) -> &'static [T]
    {
        let Self {ptr, len} = self;
        if ptr.is_null()
        {
            return &[]
        }
        unsafe {
            core::slice::from_raw_parts(ptr, len)
        }
    }
}